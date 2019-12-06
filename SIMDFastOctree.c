
//
// Created by mschnur on 11/12/19.
//

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <immintrin.h>

#include "FastOctree.h"

extern long long rdtsc();

#define MAX(x,y) (((x) > (y)) ? (x) : (y))
#define MIN(x,y) (((x) < (y)) ? (x) : (y))

#define TRUE 1
#define FALSE 0

//#define DEBUG_TRACE_FUNCTION_CALLS 1
//#define DEBUG_RAY_PARAMETER 1
//#define DEBUG_PROC_SUBTREE 1

const Node* const NO_CHILD = NULL;

const Node* const NO_CHILDREN[8] = {
        NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL
};

double CLAMPING_THRES_MIN;
double CLAMPING_THRES_MAX;
double PROB_HIT_LOG;
double PROB_MISS_LOG;
double OCC_PROB_THRES_LOG;

const double SIZE_LOOKUP_TABLE[MAX_DEPTH + 1] = {
        RESOLUTION * (double) (1u << (MAX_DEPTH - 0u)),
        RESOLUTION * (double) (1u << (MAX_DEPTH - 1u)),
        RESOLUTION * (double) (1u << (MAX_DEPTH - 2u)),
        RESOLUTION * (double) (1u << (MAX_DEPTH - 3u)),
        RESOLUTION * (double) (1u << (MAX_DEPTH - 4u)),
        RESOLUTION * (double) (1u << (MAX_DEPTH - 5u)),
        RESOLUTION * (double) (1u << (MAX_DEPTH - 6u)),
        RESOLUTION * (double) (1u << (MAX_DEPTH - 7u)),
        RESOLUTION * (double) (1u << (MAX_DEPTH - 8u)),
        RESOLUTION * (double) (1u << (MAX_DEPTH - 9u)),
        RESOLUTION * (double) (1u << (MAX_DEPTH - 10u)),
        RESOLUTION * (double) (1u << (MAX_DEPTH - 11u)),
        RESOLUTION * (double) (1u << (MAX_DEPTH - 12u)),
        RESOLUTION * (double) (1u << (MAX_DEPTH - 13u)),
        RESOLUTION * (double) (1u << (MAX_DEPTH - 14u)),
        RESOLUTION * (double) (1u << (MAX_DEPTH - 15u)),
        RESOLUTION * (double) (1u << (MAX_DEPTH - 16u))
};

static const size_t FastOctree_to_octomap_index_map[8] = {
        0, 4, 2, 6, 1, 5, 3, 7
};

//###########################################################################################
// Utility functions
//###########################################################################################

static inline int nodeHasAnyChildren(const Node* n)
{
    return (memcmp(n->children, NO_CHILDREN, 8 * sizeof(Node*)) != 0);
}

static inline int nodeIsPrunable(const Node* n)
{
    // If all of this node's children have the same log-likelihood and have no children themselves, then we can
    // prune these children, meaning we delete them.
    for (unsigned int i = 0; i < 8; ++i)
    {
        Node* child = n->children[i];
        if (child == NO_CHILD || child->logOdds != n->logOdds)
        {
            return FALSE;
        }
    }

    return TRUE;
}

static inline void deleteChild(Node* n, unsigned int childIndex)
{
    free(n->children[childIndex]);
    n->children[childIndex] = NULL;
}


static inline void deleteAllChildren(Node* n)
{
    for (unsigned int i = 0; i < 8; ++i)
    {
        deleteChild(n, i);
    }
}

static inline void createChild(Node* n, unsigned int childIndex)
{
    // In order for us to be able to free the nodes individually (which may happen during pruning) we need to
    // allocate each node separately, instead of the whole array at once. Using calloc instead of malloc
    // initializes the memory to zero, which means that the this new node's `children` array will be filled with
    // zeros (which is equivalent to filling it will NULL pointers).
    n->children[childIndex] = (Node*) calloc(1, sizeof(Node));
}

static inline void createChildIfItDoesntExist(Node* n, unsigned int childIndex)
{
    if (n->children[childIndex] == NO_CHILD)
    {
        createChild(n, childIndex);
    }
}

static inline double maxChildLogLikelihood(const Node* n)
{
    double maxLogLikelihood = -INFINITY;
    for (unsigned int i = 0; i < 8; ++i)
    {
        Node* child = n->children[i];
        if (child != NO_CHILD)
        {
            if (child->logOdds > maxLogLikelihood)
            {
                maxLogLikelihood = child->logOdds;
            }
            #ifdef DEBUG_PROC_SUBTREE
                        printf("    maxChildLogLikelihood: Child %u has log-odds of %lf\n", i, child->logOdds);
            #endif
        } else {
            #ifdef DEBUG_PROC_SUBTREE
                    printf("    maxChildLogLikelihood: No child %u\n", i);
            #endif
        }
    }

    return maxLogLikelihood;
}

static inline void expandPrunedNode(Node* n)
{
    assert(!nodeHasAnyChildren(n));

    for (unsigned int i = 0; i < 8; ++i)
    {
        createChild(n, i);

        // The children should all start with the same log odds as their parent.
        n->children[i]->logOdds = n->logOdds;
    }
}

static inline int is_less(double pointX, double pointY, double pointZ,
                   double minX, double minY, double minZ)
{
    return (minX > pointX && minY > pointY && minZ > pointZ);
}

static inline int is_greater(double pointX, double pointY, double pointZ,
                      double maxX, double maxY, double maxZ)
{
    return (pointX > maxX && pointY > maxY && pointZ > maxZ);
}

static inline int any_is_less(double pointX, double pointY, double pointZ,
                          double minX, double minY, double minZ)
{
    return (minX > pointX || minY > pointY || minZ > pointZ);
}

static inline int any_is_greater(double pointX, double pointY, double pointZ,
                             double maxX, double maxY, double maxZ)
{
    return (pointX > maxX || pointY > maxY || pointZ > maxZ);
}

static inline int is_between(double minX, double minY, double minZ,
                      double pointX, double pointY, double pointZ,
                      double maxX, double maxY, double maxZ)
{
    return (minX <= pointX && pointX <= maxX &&
            minY <= pointY && pointY <= maxY &&
            minZ <= pointZ && pointZ <= maxZ);
}

static inline int is_between_vector3d(const Vector3d* min,
                               const Vector3d* point,
                               const Vector3d* max)
{
    return (min->x <= point->x && point->x <= max->x &&
            min->y <= point->y && point->y <= max->y &&
            min->z <= point->z && point->z <= max->z);
}

//###########################################################################################
// Core algorithm functions
//###########################################################################################
const long long signbit = 0x8000000000000000ull;
const long long all_ones = -1;
const long long one = 1;
const long long node_vals[8] = {0x1,0x2,0x4,0x8,0x10,0x20,0x40,0x80};
const long long zero = 0;

static inline __m256i new_node(__m256d txmv, int x, __m256d tymv, int y, __m256d tzmv, int z)
{
    long x_shift = (63-x);
    long y_shift = (63-y);
    long z_shift = (63-z);

    __m256d tmp1v = _mm256_sub_pd(txmv, tymv);
    __m256d tmp2v = _mm256_sub_pd(txmv, tzmv);
    __m256d tmp3v = _mm256_sub_pd(tymv, tzmv);

    __m256d signbitv = _mm256_broadcast_sd((double*)&signbit);
    __m128 xshv = _mm_broadcast_ss((float*)&x_shift);
    __m128 yshv = _mm_broadcast_ss((float*)&y_shift);
    __m128 zshv = _mm_broadcast_ss((float*)&z_shift);

    //X
    //currentNode |= (unsigned long long)(*((unsigned long long*)(&tmp1)) & *((unsigned long long*)(&tmp2)) & 0x8000000000000000)>>;
    __m256d tmpv = _mm256_and_pd(tmp1v,tmp2v);
    tmpv = _mm256_and_pd(tmpv,signbitv);
    __m256i currentNodev = _mm256_srl_epi64(_mm256_castpd_si256(tmpv),_mm_castps_si128(xshv));

    
    //Y
    //currentNode |= (unsigned long long)(~*((unsigned long long*)(&tmp1)) & *((unsigned long long*)(&tmp3)) & 0x8000000000000000)>>;
    tmpv = _mm256_andnot_pd(tmp1v,tmp3v);
    tmpv = _mm256_and_pd(tmpv,signbitv);
    tmpv = _mm256_castsi256_pd(_mm256_srl_epi64(_mm256_castpd_si256(tmpv),_mm_castps_si128(yshv)));
    currentNodev = _mm256_or_si256(currentNodev, _mm256_castpd_si256(tmpv));

    
    //Z
    //currentNode |= (unsigned long long)( ~*(unsigned long long*)(&tmp3) & ~*(unsigned long long*)(&tmp2) & 0x8000000000000000 )>>;
    __m256i onesv = _mm256_castpd_si256(_mm256_broadcast_sd((double*)&all_ones));
    tmp2v = _mm256_castsi256_pd(_mm256_xor_si256(_mm256_castpd_si256(tmp2v),onesv));
    tmpv = _mm256_andnot_pd(tmp3v,tmp2v);
    tmpv = _mm256_and_pd(tmpv,signbitv);
    tmpv = _mm256_castsi256_pd(_mm256_srl_epi64(_mm256_castpd_si256(tmpv),_mm_castps_si128(zshv)));
    currentNodev = _mm256_or_si256(currentNodev, _mm256_castpd_si256(tmpv));

    return currentNodev;
}

static inline __m256i compute_valid_node(__m256d t1v,__m256d t2v,__m256d t3v, __m256d tx1v, __m256d ty1v, __m256d tz1v){
    //long shift = 63;
    //__m128 shiftv = _mm_broadcast_ss((float*)&shift);
    __m256d onesv = _mm256_broadcast_sd((double*)&all_ones);


    __m256d tmp = _mm256_and_pd(t1v,t2v);
    tx1v = _mm256_xor_pd(tx1v,onesv);
    tmp = _mm256_and_pd(tmp,t3v);
    ty1v = _mm256_xor_pd(ty1v,onesv);
    tmp = _mm256_and_pd(tmp,tx1v);
    tz1v = _mm256_xor_pd(tz1v,onesv);
    tmp = _mm256_and_pd(tmp,ty1v);
    tmp = _mm256_and_pd(tmp,tz1v);

    tmp = _mm256_castsi256_pd(_mm256_srli_epi64(_mm256_castpd_si256(tmp), 32));
    return _mm256_srai_epi32(_mm256_castpd_si256(tmp), 31);

    
    /*
    return (long long)(*((unsigned long long*)(&t1)) 
                        & *((unsigned long long*)(&t2))
                        & *((unsigned long long*)(&t3))
                        & ~*((unsigned long long*)(&tx1))
                        & ~*((unsigned long long*)(&ty1))
                        & ~*((unsigned long long*)(&tz1)))>>63;*/
}

double** proc_subtree(double* tx0, double* ty0, double* tz0,
                  double* tx1, double* ty1, double* tz1,
                  int numRays, unsigned int depth,
                  Node* n, unsigned char* a, double* endpoint) {
    if (depth == MAX_DEPTH) {
        int hits = 0;
        for(int i = 0; i < numRays; i++){
            if (!any_is_greater(endpoint[i], endpoint[i], endpoint[i], tx1[i], ty1[i], tz1[i])) {
                hits += 1;
                break; 
            }
        }

        double logLikelihoodUpdate = 0.0;
        if(hits > 0){
            logLikelihoodUpdate = PROB_HIT_LOG;
        } else {
            logLikelihoodUpdate = PROB_MISS_LOG;
        }

        // Do the update
        n->logOdds += logLikelihoodUpdate;

        // Clamp the logOdds between the min/max
        n->logOdds = fmax(CLAMPING_THRES_MIN, fmin(n->logOdds, CLAMPING_THRES_MAX));
        
        free(tx0);
        free(ty0);
        free(tz0);
        free(tx1);
        free(ty1);
        free(tz1);
        free(endpoint);
        free(a);

        return;
    }

    
    //compute t_ms
    double *txm, *tym, *tzm;
    posix_memalign((void**)&txm,64,numRays*sizeof(double));
    posix_memalign((void**)&tym,64,numRays*sizeof(double));
    posix_memalign((void**)&tzm,64,numRays*sizeof(double));

    unsigned char* nodes = (unsigned char*)calloc(numRays,sizeof(unsigned char));
    int cur_index[8] = {0};

    long long start = rdtsc();

    int diff = numRays%4;
    int i;
    for(i = 0; i < (numRays-diff); i+=4){
            /*
            int currentNode = 0;
            double t1,t2,t3,t4,t5,t6;
            int eq;
            unsigned char tmp_node = 0;
            */

            __m256d tx0v = _mm256_load_pd(tx0+i);
            __m256d ty0v = _mm256_load_pd(ty0+i);
            __m256d tz0v = _mm256_load_pd(tz0+i);
            __m256d tx1v = _mm256_load_pd(tx1+i);
            __m256d ty1v = _mm256_load_pd(ty1+i);
            __m256d tz1v = _mm256_load_pd(tz1+i);

            double half = 0.5;
            __m256d halfv = _mm256_broadcast_sd(&half);

            /*
            double txmt = 0.5 * (tx0[i] + tx1[i]);
            double tymt = 0.5 * (ty0[i] + ty1[i]);
            double tzmt = 0.5 * (tz0[i] + tz1[i]);
            */

            __m256d txmv = _mm256_add_pd(tx0v,tx1v);
            txmv = _mm256_mul_pd(txmv,halfv);
            __m256d tymv = _mm256_add_pd(ty0v,ty1v);
            tymv = _mm256_mul_pd(tymv,halfv);
            __m256d tzmv = _mm256_add_pd(tz0v,tz1v);
            tzmv = _mm256_mul_pd(tzmv,halfv);


            /*
            double tmp1 = ty0[i] - tx0[i];
            double tmp2 = tz0[i] - tx0[i];
            double tmp3 = tymt - tx0[i];
            double tmp4 = tzmt - tx0[i];
            */

            __m256d tmp1v = _mm256_sub_pd(ty0v, tx0v);
            __m256d tmp2v = _mm256_sub_pd(tz0v, tx0v);
            __m256d tmp3v = _mm256_sub_pd(tymv, tx0v);
            __m256d tmp4v = _mm256_sub_pd(tzmv, tx0v);

            /*
            currentNode |= (int)((unsigned long long)(*((unsigned long long*)(&tmp1)) & *((unsigned long long*)(&tmp2)) & *((unsigned long long*)(&tmp3)) & 0x8000000000000000)>>62);
            currentNode |= (int)((unsigned long long)(*((unsigned long long*)(&tmp1)) & *((unsigned long long*)(&tmp2)) & *((unsigned long long*)(&tmp4)) & 0x8000000000000000)>>63);
            */
           __m256i signbitv = _mm256_castpd_si256(_mm256_broadcast_sd((double*)&signbit));
           __m256i onesv = _mm256_castpd_si256(_mm256_broadcast_sd((double*)&all_ones));

            __m256i tmpv = _mm256_and_si256(_mm256_castpd_si256(tmp1v),_mm256_castpd_si256(tmp2v));
            tmpv = _mm256_and_si256(tmpv,_mm256_castpd_si256(tmp3v));
            __m256i currentNodev = _mm256_and_si256(tmpv, signbitv);
            tmpv = _mm256_and_si256(_mm256_castpd_si256(tmp1v),_mm256_castpd_si256(tmp2v));
            tmpv = _mm256_and_si256(tmpv,_mm256_castpd_si256(tmp4v));
            tmpv = _mm256_and_si256(tmpv, signbitv);
            currentNodev = _mm256_or_si256(currentNodev,tmpv);
            

            /*
            double tmp5 = tz0[i] - ty0[i];
            tmp3 = txmt - ty0[i];
            tmp4 = tzmt - ty0[i];
            */

            __m256d tmp5v = _mm256_sub_pd(tz0v, ty0v);
            tmp3v = _mm256_sub_pd(txmv, ty0v);
            tmp4v = _mm256_sub_pd(tzmv, ty0v);
            
            /*
            currentNode |= (int)((unsigned long long)(~*((unsigned long long*)(&tmp1)) & *((unsigned long long*)(&tmp5)) & *((unsigned long long*)(&tmp3)) & 0x8000000000000000)>>61);
            currentNode |= (int)((unsigned long long)(~*((unsigned long long*)(&tmp1)) & *((unsigned long long*)(&tmp5)) & *((unsigned long long*)(&tmp4)) & 0x8000000000000000)>>63);
            */
            tmpv = _mm256_andnot_si256(_mm256_castpd_si256(tmp1v),_mm256_castpd_si256(tmp5v));
            tmpv = _mm256_and_si256(tmpv,_mm256_castpd_si256(tmp3v));
            tmpv = _mm256_and_si256(tmpv, signbitv);
            currentNodev = _mm256_or_si256(currentNodev,tmpv);
            tmpv = _mm256_andnot_si256(_mm256_castpd_si256(tmp1v),_mm256_castpd_si256(tmp5v));
            tmpv = _mm256_and_si256(tmpv,_mm256_castpd_si256(tmp4v));
            tmpv = _mm256_and_si256(tmpv, signbitv);
            currentNodev = _mm256_or_si256(currentNodev,tmpv);

            /*
            tmp3 = txmt - tz0[i];
            tmp4 = tymt - tz0[i];
            */

            tmp3v = _mm256_sub_pd(txmv, tz0v);
            tmp4v = _mm256_sub_pd(tymv, tz0v);

            /*
            currentNode |= (int)((unsigned long long)(~*((unsigned long long*)(&tmp2)) & ~*((unsigned long long*)(&tmp5)) & *((unsigned long long*)(&tmp3)) & 0x8000000000000000)>>61);
            currentNode |= (int)((unsigned long long)(~*((unsigned long long*)(&tmp2)) & ~*((unsigned long long*)(&tmp5)) & *((unsigned long long*)(&tmp4)) & 0x8000000000000000)>>62);
            */
            
            tmp5v = _mm256_castsi256_pd(_mm256_xor_si256(_mm256_castpd_si256(tmp5v),onesv));
            tmpv = _mm256_andnot_si256(_mm256_castpd_si256(tmp2v),_mm256_castpd_si256(tmp5v));
            tmpv = _mm256_and_si256(tmpv,_mm256_castpd_si256(tmp3v));
            tmpv = _mm256_and_si256(tmpv, signbitv);
            currentNodev = _mm256_or_si256(currentNodev,tmpv);
            tmpv = _mm256_andnot_si256(_mm256_castpd_si256(tmp2v),_mm256_castpd_si256(tmp5v));
            tmpv = _mm256_and_si256(tmpv,_mm256_castpd_si256(tmp4v));
            tmpv = _mm256_and_si256(tmpv, signbitv);
            currentNodev = _mm256_or_si256(currentNodev,tmpv);

            onesv = _mm256_castpd_si256(_mm256_broadcast_sd((double*)&one));
            //currentNode = 1u<<currentNode;
            currentNodev = _mm256_sll_epi64(onesv, _mm256_castsi256_si128(currentNodev));


            /*
            t1 = tx0[i] - endpoint[i];
            t2 = ty0[i] - endpoint[i];
            t3 = tz0[i] - endpoint[i];
            t4 = txmt - endpoint[i];
            t5 = tymt - endpoint[i];
            t6 = tzmt - endpoint[i];
            */

            __m256d endpointv = _mm256_load_pd(endpoint+i);
            __m256d t1v = _mm256_sub_pd(tx0v, endpointv);
            __m256d t2v = _mm256_sub_pd(ty0v, endpointv);
            __m256d t3v = _mm256_sub_pd(tz0v, endpointv);
            __m256d t4v = _mm256_sub_pd(txmv, endpointv);
            __m256d t5v = _mm256_sub_pd(tymv, endpointv);
            __m256d t6v = _mm256_sub_pd(tzmv, endpointv);

            
            __m256i nextNode0, nextNode1, nextNode2, nextNode3, nextNode4, nextNode5, nextNode6;
            nextNode0 = new_node(txmv, 4, tymv, 2, tzmv, 1);
            nextNode1 = new_node(txmv, 5, tymv, 3, tz1v, 8);
            nextNode2 = new_node(txmv, 6, ty1v, 8, tzmv, 3);       
            nextNode3 = new_node(txmv, 7, ty1v, 8, tz1v, 8);
            nextNode4 = new_node(tx1v, 8, tymv, 6, tzmv, 5);
            nextNode5 = new_node(tx1v, 8, tymv, 7, tz1v, 8);
            nextNode6 = new_node(tx1v, 8, ty1v, 8, tzmv, 7);

            /*
            eq = ~(((1<<0) - currentNode)>>31);
            tmp_node = (unsigned char)(currentNode & valid_node & eq);
            cur_index[a[i]] += (1 & valid_node & eq);
            currentNode = (nextNode0 & eq) | (currentNode & ~eq);
            */
            __m256i tmp_av = _mm256_castpd_si256(_mm256_broadcast_sd((double*)&zero));
            onesv = _mm256_castpd_si256(_mm256_broadcast_sd((double*)&all_ones));
            __m256i neqv = _mm256_castpd_si256(_mm256_broadcast_sd((double*)node_vals));
            neqv = _mm256_sub_epi64(neqv,currentNodev);
            neqv = _mm256_srli_epi64(neqv, 32);
            neqv = _mm256_srai_epi32(neqv, 31);
            __m256i eqv = _mm256_xor_si256(neqv,onesv);
            __m256i valid_nodev = compute_valid_node(t1v,t2v,t3v,txmv,tymv,tzmv);
            __m256i tmp_nodev = _mm256_and_si256(valid_nodev,eqv);
            tmp_nodev = _mm256_and_si256(tmp_nodev, currentNodev);
            tmp_av = _mm256_or_si256(tmp_av,tmp_nodev);
            tmpv = _mm256_and_si256(currentNodev,neqv);
            currentNodev = _mm256_and_si256(nextNode0,eqv);
            currentNodev = _mm256_or_si256(currentNodev,tmpv);

            /*
            eq = ~(((1<<1) - currentNode)>>31);
            valid_node = compute_valid_node(t1v,t2v,t6v,txmv,tymv,tz1v);
            tmp_node |= (unsigned char)(currentNode & valid_node & eq);
            cur_index[1u^a[i]] += (1 & valid_node & eq);
            currentNode = (nextNode1 & eq) | (currentNode & ~eq);
            */
            neqv = _mm256_castpd_si256(_mm256_broadcast_sd((double*)(node_vals+1)));
            neqv = _mm256_sub_epi64(neqv,currentNodev);
            neqv = _mm256_srli_epi64(neqv, 32);
            neqv = _mm256_srai_epi32(neqv, 31);
            eqv = _mm256_xor_si256(neqv,onesv);
            valid_nodev = compute_valid_node(t1v,t2v,t6v,txmv,tymv,tz1v);
            tmp_nodev = _mm256_and_si256(valid_nodev,eqv);
            tmp_nodev = _mm256_and_si256(tmp_nodev, currentNodev);
            tmp_av = _mm256_or_si256(tmp_av,tmp_nodev);
            tmpv = _mm256_and_si256(currentNodev,neqv);
            currentNodev = _mm256_and_si256(nextNode1,eqv);
            currentNodev = _mm256_or_si256(currentNodev,tmpv);

            /*
            eq = ~(((1<<2) - currentNode)>>31);
            valid_node = compute_valid_node(t1v,t5v,t3v,txmv,ty1v,tzmv);
            tmp_node |= (unsigned char)(currentNode & valid_node & eq);
            cur_index[2u^a[i]] += (1 & valid_node & eq);
            currentNode = (nextNode2 & eq) | (currentNode & ~eq);
            */
            neqv = _mm256_castpd_si256(_mm256_broadcast_sd((double*)(node_vals+2)));
            neqv = _mm256_sub_epi64(neqv,currentNodev);
            neqv = _mm256_srli_epi64(neqv, 32);
            neqv = _mm256_srai_epi32(neqv, 31);
            eqv = _mm256_xor_si256(neqv,onesv);
            valid_nodev = compute_valid_node(t1v,t5v,t3v,txmv,ty1v,tzmv);
            tmp_nodev = _mm256_and_si256(valid_nodev,eqv);
            tmp_nodev = _mm256_and_si256(tmp_nodev, currentNodev);
            tmp_av = _mm256_or_si256(tmp_av,tmp_nodev);
            tmpv = _mm256_and_si256(currentNodev,neqv);
            currentNodev = _mm256_and_si256(nextNode2,eqv);
            currentNodev = _mm256_or_si256(currentNodev,tmpv);

            /*
            eq = ~(((1<<3) - currentNode)>>31);
            valid_node = compute_valid_node(t1v,t5v,t6v,txm,ty1v,tz1v);
            tmp_node |= (unsigned char)(currentNode & valid_node & eq);
            cur_index[3u^a[i]] += (1 & valid_node & eq);
            currentNode = (nextNode3  & eq) | (currentNode & ~eq);
            */
            neqv = _mm256_castpd_si256(_mm256_broadcast_sd((double*)(node_vals+3)));
            neqv = _mm256_sub_epi64(neqv,currentNodev);
            neqv = _mm256_srli_epi64(neqv, 32);
            neqv = _mm256_srai_epi32(neqv, 31);
            eqv = _mm256_xor_si256(neqv,onesv);
            valid_nodev = compute_valid_node(t1v,t5v,t6v,txmv,ty1v,tz1v);
            tmp_nodev = _mm256_and_si256(valid_nodev,eqv);
            tmp_nodev = _mm256_and_si256(tmp_nodev, currentNodev);
            tmp_av = _mm256_or_si256(tmp_av,tmp_nodev);
            tmpv = _mm256_and_si256(currentNodev,neqv);
            currentNodev = _mm256_and_si256(nextNode3,eqv);
            currentNodev = _mm256_or_si256(currentNodev,tmpv);

            /*
            eq = ~(((1<<4) - currentNode)>>31);
            valid_node = compute_valid_node(t4v,t2v,t3v,tx1v,tymv,tzmv);
            tmp_node |= (unsigned char)(currentNode & valid_node & eq);
            cur_index[4u^a[i]] += (1 & valid_node & eq);
            currentNode = (nextNode4 & eq) | (currentNode & ~eq);
            */
            neqv = _mm256_castpd_si256(_mm256_broadcast_sd((double*)(node_vals+4)));
            neqv = _mm256_sub_epi64(neqv,currentNodev);
            neqv = _mm256_srli_epi64(neqv, 32);
            neqv = _mm256_srai_epi32(neqv, 31);
            eqv = _mm256_xor_si256(neqv,onesv);
            valid_nodev = compute_valid_node(t4v,t2v,t3v,tx1v,tymv,tzmv);
            tmp_nodev = _mm256_and_si256(valid_nodev,eqv);
            tmp_nodev = _mm256_and_si256(tmp_nodev, currentNodev);
            tmp_av = _mm256_or_si256(tmp_av,tmp_nodev);
            tmpv = _mm256_and_si256(currentNodev,neqv);
            currentNodev = _mm256_and_si256(nextNode4,eqv);
            currentNodev = _mm256_or_si256(currentNodev,tmpv);

            /*
            eq = ~(((1<<5) - currentNode)>>31);
            valid_node = compute_valid_node(t4v,t2v,t6v,tx1v,tymv,tz1v);
            tmp_node |= (unsigned char)(currentNode & valid_node & eq);
            cur_index[5u^a[i]] += (1 & valid_node & eq);
            currentNode = (nextNode5 & eq) | (currentNode & ~eq);
            */
            neqv = _mm256_castpd_si256(_mm256_broadcast_sd((double*)(node_vals+5)));
            neqv = _mm256_sub_epi64(neqv,currentNodev);
            neqv = _mm256_srli_epi64(neqv, 32);
            neqv = _mm256_srai_epi32(neqv, 31);
            eqv = _mm256_xor_si256(neqv,onesv);
            valid_nodev = compute_valid_node(t4v,t2v,t6v,tx1v,tymv,tz1v);
            tmp_nodev = _mm256_and_si256(valid_nodev,eqv);
            tmp_nodev = _mm256_and_si256(tmp_nodev, currentNodev);
            tmp_av = _mm256_or_si256(tmp_av,tmp_nodev);
            tmpv = _mm256_and_si256(currentNodev,neqv);
            currentNodev = _mm256_and_si256(nextNode5,eqv);
            currentNodev = _mm256_or_si256(currentNodev,tmpv);

           /*
            eq = ~(((1<<6) - currentNode)>>31);
            valid_node = compute_valid_node(t4v,t5v,t3v,tx1v,ty1v,tzmv);
            tmp_node |= (unsigned char)(currentNode & valid_node & eq);
            cur_index[6u^a[i]] += (1 & valid_node & eq);
            currentNode = (nextNode6 & eq) | (currentNode & ~eq);
            */
            neqv = _mm256_castpd_si256(_mm256_broadcast_sd((double*)(node_vals+6)));
            neqv = _mm256_sub_epi64(neqv,currentNodev);
            neqv = _mm256_srli_epi64(neqv, 32);
            neqv = _mm256_srai_epi32(neqv, 31);
            eqv = _mm256_xor_si256(neqv,onesv);
            valid_nodev = compute_valid_node(t4v,t5v,t3v,tx1v,ty1v,tzmv);
            tmp_nodev = _mm256_and_si256(valid_nodev,eqv);
            tmp_nodev = _mm256_and_si256(tmp_nodev, currentNodev);
            tmp_av = _mm256_or_si256(tmp_av,tmp_nodev);
            tmpv = _mm256_and_si256(currentNodev,neqv);
            currentNodev = _mm256_and_si256(nextNode6,eqv);
            currentNodev = _mm256_or_si256(currentNodev,tmpv);

           /*
            eq = ~(((1<<7) - currentNode)>>31);
            valid_node = compute_valid_node(t4v,t5v,t6v,tx1v,ty1v,tz1v);
            tmp_node |= (unsigned char)(currentNode & valid_node & eq);
            nodes[i] |= tmp_node;
            cur_index[7u^a[i]] += (1 & valid_node & eq);
            */
            neqv = _mm256_castpd_si256(_mm256_broadcast_sd((double*)(node_vals+7)));
            neqv = _mm256_sub_epi64(neqv,currentNodev);
            neqv = _mm256_srli_epi64(neqv, 32);
            neqv = _mm256_srai_epi32(neqv, 31);
            eqv = _mm256_xor_si256(neqv,onesv);
            valid_nodev = compute_valid_node(t4v,t5v,t6v,tx1v,ty1v,tz1v);
            tmp_nodev = _mm256_and_si256(valid_nodev,eqv);
            tmp_nodev = _mm256_and_si256(tmp_nodev, currentNodev);
            tmp_av = _mm256_or_si256(tmp_av,tmp_nodev);
            
            /*
            txm[i] = txmt;
            tym[i] = tymt;
            tzm[i] = tzmt;
            */
            _mm256_store_pd((txm+i), txmv);
            _mm256_store_pd((tym+i), tymv);
            _mm256_store_pd((tzm+i), tzmv);
    }

    
    long long end = rdtsc();
    printf("proc_subtree cycles: %d\n", (end-start));
    printf("%d\n",i);

    //so the loop above doesn't get optimized away
    double** tms = (double**)malloc(3*sizeof(double**));

    tms[0] = txm;
    tms[1] = tym;
    tms[2] = tzm;

    return tms;

/*
    double* new_tx0[8] = {0};
    double* new_ty0[8] = {0};
    double* new_tz0[8] = {0};
    double* new_tx1[8] = {0};
    double* new_ty1[8] = {0};
    double* new_tz1[8] = {0};
    double* new_endpoints[8] = {0};
    unsigned char* new_a[8] = {0};

    for(int i = 0; i<8; i++){
        if(cur_index[i] > 0){
            new_tx0[i] = (double*)malloc(cur_index[i]*sizeof(double));
            new_ty0[i] = (double*)malloc(cur_index[i]*sizeof(double));
            new_tz0[i] = (double*)malloc(cur_index[i]*sizeof(double));
            new_tx1[i] = (double*)malloc(cur_index[i]*sizeof(double));
            new_ty1[i] = (double*)malloc(cur_index[i]*sizeof(double));
            new_tz1[i] = (double*)malloc(cur_index[i]*sizeof(double));
            new_endpoints[i] = (double*)malloc(cur_index[i]*sizeof(double));
            new_a[i] = (unsigned char*)calloc(cur_index[i], sizeof(unsigned char));
        }
    }

    int update_index[8] = {0};
    for(int i = 0; i < numRays; i++){
        if(nodes[i] & (1u<<0)){
            new_tx0[a[i]][update_index[a[i]]] = tx0[i];
            new_ty0[a[i]][update_index[a[i]]] = ty0[i];
            new_tz0[a[i]][update_index[a[i]]] = tz0[i];
            new_tx1[a[i]][update_index[a[i]]] = txm[i];
            new_ty1[a[i]][update_index[a[i]]] = tym[i];
            new_tz1[a[i]][update_index[a[i]]] = tzm[i];
            new_endpoints[a[i]][update_index[a[i]]] = endpoint[i];
            new_a[a[i]][update_index[a[i]]] = a[i];
            update_index[a[i]]++;
        }
        if(nodes[i] & (1u<<1)){
            new_tx0[1u^a[i]][update_index[1u^a[i]]] = tx0[i];
            new_ty0[1u^a[i]][update_index[1u^a[i]]] = ty0[i];
            new_tz0[1u^a[i]][update_index[1u^a[i]]] = tzm[i];
            new_tx1[1u^a[i]][update_index[1u^a[i]]] = txm[i];
            new_ty1[1u^a[i]][update_index[1u^a[i]]] = tym[i];
            new_tz1[1u^a[i]][update_index[1u^a[i]]] = tz1[i];
            new_endpoints[1u^a[i]][update_index[1u^a[i]]] = endpoint[i];
            new_a[1u^a[i]][update_index[1u^a[i]]] = a[i];
            update_index[1u^a[i]]++;
        }
        if(nodes[i] & (1u<<2)){
            new_tx0[2u^a[i]][update_index[2u^a[i]]] = tx0[i];
            new_ty0[2u^a[i]][update_index[2u^a[i]]] = tym[i];
            new_tz0[2u^a[i]][update_index[2u^a[i]]] = tz0[i];
            new_tx1[2u^a[i]][update_index[2u^a[i]]] = txm[i];
            new_ty1[2u^a[i]][update_index[2u^a[i]]] = ty1[i];
            new_tz1[2u^a[i]][update_index[2u^a[i]]] = tzm[i];
            new_endpoints[2u^a[i]][update_index[2u^a[i]]] = endpoint[i];
            new_a[2u^a[i]][update_index[2u^a[i]]] = a[i];
            update_index[2u^a[i]]++;
        }
        if(nodes[i] & (1u<<3)){
            new_tx0[3u^a[i]][update_index[3u^a[i]]] = tx0[i];
            new_ty0[3u^a[i]][update_index[3u^a[i]]] = tym[i];
            new_tz0[3u^a[i]][update_index[3u^a[i]]] = tzm[i];
            new_tx1[3u^a[i]][update_index[3u^a[i]]] = txm[i];
            new_ty1[3u^a[i]][update_index[3u^a[i]]] = ty1[i];
            new_tz1[3u^a[i]][update_index[3u^a[i]]] = tz1[i];
            new_endpoints[3u^a[i]][update_index[3u^a[i]]] = endpoint[i];
            new_a[3u^a[i]][update_index[3u^a[i]]] = a[i];
            update_index[3u^a[i]]++;
        }
        if(nodes[i] & (1u<<4)){
            new_tx0[4u^a[i]][update_index[4u^a[i]]] = txm[i];
            new_ty0[4u^a[i]][update_index[4u^a[i]]] = ty0[i];
            new_tz0[4u^a[i]][update_index[4u^a[i]]] = tz0[i];
            new_tx1[4u^a[i]][update_index[4u^a[i]]] = tx1[i];
            new_ty1[4u^a[i]][update_index[4u^a[i]]] = tym[i];
            new_tz1[4u^a[i]][update_index[4u^a[i]]] = tzm[i];
            new_endpoints[4u^a[i]][update_index[4u^a[i]]] = endpoint[i];
            new_a[4u^a[i]][update_index[4u^a[i]]] = a[i];
            update_index[4u^a[i]]++;
        }
        if(nodes[i] & (1u<<5)){
            new_tx0[5u^a[i]][update_index[5u^a[i]]] = txm[i];
            new_ty0[5u^a[i]][update_index[5u^a[i]]] = ty0[i];
            new_tz0[5u^a[i]][update_index[5u^a[i]]] = tzm[i];
            new_tx1[5u^a[i]][update_index[5u^a[i]]] = tx1[i];
            new_ty1[5u^a[i]][update_index[5u^a[i]]] = tym[i];
            new_tz1[5u^a[i]][update_index[5u^a[i]]] = tz1[i];
            new_endpoints[5u^a[i]][update_index[5u^a[i]]] = endpoint[i];
            new_a[5u^a[i]][update_index[5u^a[i]]] = a[i];
            update_index[5u^a[i]]++;
        }
        if(nodes[i] & (1u<<6)){
            new_tx0[6u^a[i]][update_index[6u^a[i]]] = txm[i];
            new_ty0[6u^a[i]][update_index[6u^a[i]]] = tym[i];
            new_tz0[6u^a[i]][update_index[6u^a[i]]] = tz0[i];
            new_tx1[6u^a[i]][update_index[6u^a[i]]] = tx1[i];
            new_ty1[6u^a[i]][update_index[6u^a[i]]] = ty1[i];
            new_tz1[6u^a[i]][update_index[6u^a[i]]] = tzm[i];
            new_endpoints[6u^a[i]][update_index[6u^a[i]]] = endpoint[i];
            new_a[6u^a[i]][update_index[6u^a[i]]] = a[i];
            update_index[6u^a[i]]++;
        }
        if(nodes[i] & (1u<<7)){
            new_tx0[7u^a[i]][update_index[7u^a[i]]] = txm[i];
            new_ty0[7u^a[i]][update_index[7u^a[i]]] = tym[i];
            new_tz0[7u^a[i]][update_index[7u^a[i]]] = tzm[i];
            new_tx1[7u^a[i]][update_index[7u^a[i]]] = tx1[i];
            new_ty1[7u^a[i]][update_index[7u^a[i]]] = ty1[i];
            new_tz1[7u^a[i]][update_index[7u^a[i]]] = tz1[i];
            new_endpoints[7u^a[i]][update_index[7u^a[i]]] = endpoint[i];
            new_a[7u^a[i]][update_index[7u^a[i]]] = a[i];
            update_index[7u^a[i]]++;
        }         
    }

    free(tx0);
    free(ty0);
    free(tz0);
    free(tx1);
    free(ty1);
    free(tz1);
    free(endpoint);
    free(a);   

    free(nodes);

    /*for(int node = 0; node < 8; node++) {
        if(cur_index[node] > 0){
            //assumes a is already handled by now
            createChildIfItDoesntExist(n, node);
            proc_subtree(new_tx0[node], new_ty0[node], new_tz0[node],
                         new_tx1[node], new_ty1[node], new_tz1[node],
                         cur_index[node], depth + 1, 
                         n->children[node], new_a[node], new_endpoints[node]);
        }
    }
        n->logOdds = maxChildLogLikelihood(n);*/
}


void ray_parameter(Octree* tree, Ray* rays, int numRays) {
    int createdRoot = FALSE;
    if (tree->root == NULL) {
        // Using calloc instead of malloc initializes the memory to zero, which means that the the root's `children`
        // array will be filled with zeros (which is equivalent to filling it will NULL pointers). It also sets the
        // root's log odds to 0, which we want as our initial value.
        tree->root = (Node *) calloc(1, sizeof(Node));
        createdRoot = TRUE;
    }

    double *tx0, *ty0, *tz0, *tx1, *ty1, *tz1, *endpoints;
    posix_memalign((void**)&tx0,64,numRays*sizeof(double));
    posix_memalign((void**)&ty0,64,numRays*sizeof(double));
    posix_memalign((void**)&tz0,64,numRays*sizeof(double));
    posix_memalign((void**)&tx1,64,numRays*sizeof(double));
    posix_memalign((void**)&ty1,64,numRays*sizeof(double));
    posix_memalign((void**)&tz1,64,numRays*sizeof(double));
    posix_memalign((void**)&endpoints,64,numRays*sizeof(double));
    unsigned char* a = (unsigned char*)calloc(numRays, sizeof(unsigned char));

    for(int i = 0; i < numRays; i++){
        long long reflected, cur, tmp, sign_bit_if_negative;
        long long reflected1, cur1, tmp1, sign_bit_if_negative1;
        long long reflected2, cur2, tmp2, sign_bit_if_negative2;
        Ray* r = (rays+i);

        // sign_bit_if_neg = directionX & 0x8000000000000000;
        // originX = originX XOR sign_bit_if_neg;
        // directionX = directionX & 0x7FFFFFFFFFFFFFFF;
        // sign_bit_if_neg = sign_bit_if_neg >> 63; // right shift sign-extend
        // partial_a = a_number & sign_bit_if_neg;
        // local_a = local_a | partial_a;

        cur = *((unsigned long long*)(&(r->direction.x)));
        cur1 = *((unsigned long long*)(&(r->direction.y)));
        cur2 = *((unsigned long long*)(&(r->direction.z)));
        sign_bit_if_negative = cur & 0x8000000000000000ull;
        sign_bit_if_negative1 = cur1 & 0x8000000000000000ull;
        sign_bit_if_negative2 = cur2 & 0x8000000000000000ull;
        tmp = *((unsigned long long*)(&(r->origin.x)));
        tmp1 = *((unsigned long long*)(&(r->origin.y)));
        tmp2 = *((unsigned long long*)(&(r->origin.z)));
        tmp ^= sign_bit_if_negative;
        tmp1 ^= sign_bit_if_negative1;
        tmp2 ^= sign_bit_if_negative2;
        r->origin.x = *((double*)(&tmp));
        r->origin.y = *((double*)(&tmp1));
        r->origin.z = *((double*)(&tmp2));
        tmp = cur ^ sign_bit_if_negative;
        tmp1 = cur1 ^ sign_bit_if_negative1;
        tmp2 = cur2 ^ sign_bit_if_negative2;
        r->direction.x = *((double*)(&tmp));
        r->direction.y = *((double*)(&tmp1));
        r->direction.z = *((double*)(&tmp2));
        reflected = *((long long*)(&sign_bit_if_negative))>>63;
        reflected1 = *((long long*)(&sign_bit_if_negative1))>>63;
        reflected2 = *((long long*)(&sign_bit_if_negative2))>>63;
        a[i] |= (4u & reflected);
        a[i] |= (2u & reflected1);
        a[i] |= (1u & reflected2);

        // Improve IEEE double stability
        double rdxInverse = 1.0 / r->direction.x;
        double rdyInverse = 1.0 / r->direction.y;
        double rdzInverse = 1.0 / r->direction.z;

        tx0[i] = (tree->min.x - r->origin.x) * rdxInverse;
        tx1[i] = (tree->max.x - r->origin.x) * rdxInverse;
        ty0[i] = (tree->min.y - r->origin.y) * rdyInverse;
        ty1[i] = (tree->max.y - r->origin.y) * rdyInverse;
        tz0[i] = (tree->min.z - r->origin.z) * rdzInverse;
        tz1[i] = (tree->max.z - r->origin.z) * rdzInverse;
        endpoints[i] = r->t_end;
    }

    free(rays);

    // for now assume our point cloud origin and all points exist within the actree bounds
    double** mids = proc_subtree(tx0, ty0, tz0, tx1, ty1, tz1, numRays, 0, tree->root, a, endpoints);
        
}

void insertPointCloud(Octree* tree, Vector3d* points, size_t numPoints, Vector3d* sensorOrigin)
{
    static int notInitialized = TRUE;
    if (notInitialized) {
        CLAMPING_THRES_MIN = logodds(0.1192);
        CLAMPING_THRES_MAX = logodds(0.971);
        PROB_HIT_LOG = logodds(0.7);
        PROB_MISS_LOG = logodds(0.4);
        OCC_PROB_THRES_LOG = logodds(0.5);
        notInitialized = FALSE;
    }

    Ray* rays = (Ray*)calloc(numPoints, sizeof(Ray));

    // for each point, create ray and call ray_parameter. After all points done, prune the tree
    for (size_t i = 0; i < numPoints; ++i)
    {
        initRay((rays+i),
                sensorOrigin->x, sensorOrigin->y, sensorOrigin->z,
                points[i].x, points[i].y, points[i].z);
    }

    ray_parameter(tree, rays, numPoints);
}