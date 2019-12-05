
//
// Created by mschnur on 11/12/19.
//

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "FastOctree.h"

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

static inline int new_node(double txm, int x, double tym, int y, double tzm, int z)
{
    unsigned long long currentNode = 0;

    double tmp1 = txm - tym;
    double tmp2 = txm - tzm;
    double tmp3 = tym - tzm;

    //X
    currentNode |= (unsigned long long)(*((unsigned long long*)(&tmp1)) & *((unsigned long long*)(&tmp2)) & 0x8000000000000000)>>(63-x);
    
    //Y
    currentNode |= (unsigned long long)(~*((unsigned long long*)(&tmp1)) & *((unsigned long long*)(&tmp3)) & 0x8000000000000000)>>(63-y);
    
    //Z
    currentNode |= (unsigned long long)( ~*(unsigned long long*)(&tmp3) & ~*(unsigned long long*)(&tmp2) & 0x8000000000000000 )>>(63-z);

    return (int)currentNode;
}

static inline long long compute_valid_node(double t1,double t2,double t3, double tx1, double ty1, double tz1){
    return (long long)(*((unsigned long long*)(&t1)) 
                        & *((unsigned long long*)(&t2))
                        & *((unsigned long long*)(&t3))
                        & ~*((unsigned long long*)(&tx1))
                        & ~*((unsigned long long*)(&ty1))
                        & ~*((unsigned long long*)(&tz1))
                        & 0x8000000000000000)>>63;
}

int printed = 0;

void proc_subtree(double* tx0, double* ty0, double* tz0,
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
    double* txm = (double*)calloc(numRays, sizeof(double));
    double* tym = (double*)calloc(numRays, sizeof(double));
    double* tzm = (double*)calloc(numRays, sizeof(double));

    unsigned char* nodes = (unsigned char*)calloc(numRays,sizeof(unsigned char));
    int cur_index[8] = {0};
    for(int i = 0; i < numRays; i++){
        int currentNode = 0;
        double t1,t2,t3,t4,t5,t6;
        int eq;
        long long valid_node;
        unsigned char tmp_node = 0;
        

        #ifdef DEBUG_PROC_SUBTREE
            printf("txm = %lf, tym = %lf, tzm = %lf\n", txm, tym, tzm);
        #endif

        txm[i] = 0.5 * (tx0[i] + tx1[i]);
        tym[i] = 0.5 * (ty0[i] + ty1[i]);
        tzm[i] = 0.5 * (tz0[i] + tz1[i]);

        double tmp1 = ty0[i] - tx0[i];
        double tmp2 = tz0[i] - tx0[i];
        double tmp3 = tym[i] - tx0[i];
        double tmp4 = tzm[i] - tx0[i];

        double tmp5 = tz0[i] - ty0[i];
        double tmp6 = txm[i] - ty0[i];
        double tmp7 = tzm[i] - ty0[i];

        double tmp8 = txm[i] - tz0[i];
        double tmp9 = tym[i] - tz0[i];


        currentNode |= (int)((unsigned long long)(*((unsigned long long*)(&tmp1)) & *((unsigned long long*)(&tmp2)) & *((unsigned long long*)(&tmp3)) & 0x8000000000000000)>>62);
        currentNode |= (int)((unsigned long long)(*((unsigned long long*)(&tmp1)) & *((unsigned long long*)(&tmp2)) & *((unsigned long long*)(&tmp4)) & 0x8000000000000000)>>63);

        currentNode |= (int)((unsigned long long)(~*((unsigned long long*)(&tmp1)) & *((unsigned long long*)(&tmp5)) & *((unsigned long long*)(&tmp6)) & 0x8000000000000000)>>61);
        currentNode |= (int)((unsigned long long)(~*((unsigned long long*)(&tmp1)) & *((unsigned long long*)(&tmp5)) & *((unsigned long long*)(&tmp7)) & 0x8000000000000000)>>63);

        currentNode |= (int)((unsigned long long)(~*((unsigned long long*)(&tmp2)) & ~*((unsigned long long*)(&tmp5)) & *((unsigned long long*)(&tmp8)) & 0x8000000000000000)>>61);
        currentNode |= (int)((unsigned long long)(~*((unsigned long long*)(&tmp2)) & ~*((unsigned long long*)(&tmp5)) & *((unsigned long long*)(&tmp9)) & 0x8000000000000000)>>62);


        #ifdef DEBUG_PROC_SUBTREE
            printf("depth=%u: first_node: %d\n", depth, currentNode);
        #endif

        currentNode = 1u<<currentNode;

        t1 = tx0[i] - endpoint[i];
        t2 = ty0[i] - endpoint[i];
        t3 = tz0[i] - endpoint[i];
        t4 = txm[i] - endpoint[i];
        t5 = tym[i] - endpoint[i];
        t6 = tzm[i] - endpoint[i];


        eq = ~(((1<<0) - currentNode)>>31);
        //this should compute if the currentNode is valid.
        valid_node = compute_valid_node(t1,t2,t3,txm[i],tym[i],tzm[i]);
        tmp_node = (unsigned char)(currentNode & valid_node & eq);
        nodes[i] |= tmp_node;
        cur_index[a[i]] += (1 & valid_node & eq);
        currentNode = (new_node(txm[i], 4, tym[i], 2, tzm[i], 1) & eq) | (currentNode & ~eq);

        eq = ~(((1<<1) - currentNode)>>31);
        valid_node = compute_valid_node(t1,t2,t6,txm[i],tym[i],tz1[i]);
        tmp_node = (unsigned char)(currentNode & valid_node & eq);
        nodes[i] |= tmp_node;
        cur_index[1u^a[i]] += (1 & valid_node & eq);
        currentNode = (new_node(txm[i], 5, tym[i], 3, tz1[i], 8) & eq) | (currentNode & ~eq);

        eq = ~(((1<<2) - currentNode)>>31);
        valid_node = compute_valid_node(t1,t5,t3,txm[i],ty1[i],tzm[i]);
        tmp_node = (unsigned char)(currentNode & valid_node & eq);
        nodes[i] |= tmp_node;
        cur_index[2u^a[i]] += (1 & valid_node & eq);
        currentNode = (new_node(txm[i], 6, ty1[i], 8, tzm[i], 3) & eq) | (currentNode & ~eq);

        eq = ~(((1<<3) - currentNode)>>31);
        valid_node = compute_valid_node(t1,t5,t6,txm[i],ty1[i],tz1[i]);
        tmp_node = (unsigned char)(currentNode & valid_node & eq);
        nodes[i] |= tmp_node;
        cur_index[3u^a[i]] += (1 & valid_node & eq);
        currentNode = (new_node(txm[i], 7, ty1[i], 8, tz1[i], 8) & eq) | (currentNode & ~eq);

        eq = ~(((1<<4) - currentNode)>>31);
        valid_node = compute_valid_node(t4,t2,t3,tx1[i],tym[i],tzm[i]);
        tmp_node = (unsigned char)(currentNode & valid_node & eq);
        nodes[i] |= tmp_node;
        cur_index[4u^a[i]] += (1 & valid_node & eq);
        currentNode = (new_node(tx1[i], 8, tym[i], 6, tzm[i], 5) & eq) | (currentNode & ~eq);

        eq = ~(((1<<5) - currentNode)>>31);
        valid_node = compute_valid_node(t4,t2,t6,tx1[i],tym[i],tz1[i]);
        tmp_node = (unsigned char)(currentNode & valid_node & eq);
        nodes[i] |= tmp_node;
        cur_index[5u^a[i]] += (1 & valid_node & eq);
        currentNode = (new_node(tx1[i], 8, tym[i], 7, tz1[i], 8) & eq) | (currentNode & ~eq);

        eq = ~(((1<<6) - currentNode)>>31);
        valid_node = compute_valid_node(t4,t5,t3,tx1[i],ty1[i],tzm[i]);
        tmp_node = (unsigned char)(currentNode & valid_node & eq);
        nodes[i] |= tmp_node;
        cur_index[6u^a[i]] += (1 & valid_node & eq);
        currentNode = (new_node(tx1[i], 8, ty1[i], 8, tzm[i], 7) & eq) | (currentNode & ~eq);

        eq = ~(((1<<7) - currentNode)>>31);
        valid_node = compute_valid_node(t4,t5,t6,tx1[i],ty1[i],tz1[i]);
        tmp_node = (unsigned char)(currentNode & valid_node & eq);
        nodes[i] |= tmp_node;
        cur_index[7u^a[i]] += (1 & valid_node & eq);
    }


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
            new_tx0[i] = (double*)calloc(cur_index[i], sizeof(double));
            new_ty0[i] = (double*)calloc(cur_index[i], sizeof(double));
            new_tz0[i] = (double*)calloc(cur_index[i], sizeof(double));
            new_tx1[i] = (double*)calloc(cur_index[i], sizeof(double));
            new_ty1[i] = (double*)calloc(cur_index[i], sizeof(double));
            new_tz1[i] = (double*)calloc(cur_index[i], sizeof(double));
            new_endpoints[i] = (double*)calloc(cur_index[i], sizeof(double));
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


    for(int node = 0; node < 8; node++) {
        if(cur_index[node] > 0){
            //assumes a is already handled by now
            createChildIfItDoesntExist(n, node);
            proc_subtree(new_tx0[node], new_ty0[node], new_tz0[node],
                         new_tx1[node], new_ty1[node], new_tz1[node],
                         cur_index[node], depth + 1, 
                         n->children[node], new_a[node], new_endpoints[node]);
        }
    }
        n->logOdds = maxChildLogLikelihood(n);
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

    double* tx0 = (double*)calloc(numRays, sizeof(double));
    double* ty0 = (double*)calloc(numRays, sizeof(double));
    double* tz0 = (double*)calloc(numRays, sizeof(double));
    double* tx1 = (double*)calloc(numRays, sizeof(double));
    double* ty1 = (double*)calloc(numRays, sizeof(double));
    double* tz1 = (double*)calloc(numRays, sizeof(double));
    double* endpoints = (double*)calloc(numRays, sizeof(double));
    unsigned char* a = (unsigned char*)calloc(numRays, sizeof(unsigned char));

    for(int i = 0; i < numRays; i++){
        long long reflected, neg, cur, tmp, sign_bit_if_negative;
        Ray* r = (rays+i);

        // sign_bit_if_neg = directionX & 0x8000000000000000;
        // originX = originX XOR sign_bit_if_neg;
        // directionX = directionX & 0x7FFFFFFFFFFFFFFF;
        // sign_bit_if_neg = sign_bit_if_neg >> 63; // right shift sign-extend
        // partial_a = a_number & sign_bit_if_neg;
        // local_a = local_a | partial_a;

        cur = *((unsigned long long*)(&(r->direction.x)));
        sign_bit_if_negative = cur & 0x8000000000000000ull;
        tmp = *((unsigned long long*)(&(r->origin.x)));
        tmp ^= sign_bit_if_negative;
        r->origin.x = *((double*)(&tmp));
        tmp = cur & 0x7fffffffffffffffull;
        r->direction.x = *((double*)(&tmp));
        reflected = *((long long*)(&sign_bit_if_negative))>>63;
        a[i] |= (4u & reflected);

        cur = *((unsigned long long*)(&(r->direction.y)));
        sign_bit_if_negative = cur & 0x8000000000000000ull;
        tmp = *((unsigned long long*)(&(r->origin.y)));
        tmp ^= sign_bit_if_negative;
        r->origin.y = *((double*)(&tmp));
        tmp = cur & 0x7fffffffffffffffull;
        r->direction.y = *((double*)(&tmp));
        reflected = *((long long*)(&sign_bit_if_negative))>>63;
        a[i] |= (2u & reflected);

        cur = *((unsigned long long*)(&(r->direction.z)));
        sign_bit_if_negative = cur & 0x8000000000000000ull;
        tmp = *((unsigned long long*)(&(r->origin.z)));
        tmp ^= sign_bit_if_negative;
        r->origin.z = *((double*)(&tmp));
        tmp = cur & 0x7fffffffffffffffull;
        r->direction.z = *((double*)(&tmp));
        reflected = *((long long*)(&sign_bit_if_negative))>>63;
        a[i] |= (1u & reflected);

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
    proc_subtree(tx0, ty0, tz0, tx1, ty1, tz1, numRays, 0, tree->root, a, endpoints);
        
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