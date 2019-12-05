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

static inline int createChildIfItDoesntExist(Node* n, unsigned int childIndex)
{
    if (n->children[childIndex] == NO_CHILD)
    {
        createChild(n, childIndex);
        return TRUE;
    }

    return FALSE;
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

int first_node(double tx0, double ty0, double tz0, double txm, double tym, double tzm)
{
    unsigned char answer = 0;	// initialize to 00000000
    // select the entry plane and set bits
    if(tx0 > ty0)
    {
        if(tx0 > tz0)
        {   // PLANE YZ
            if(tym < tx0) answer|=2u;	// set bit at position 1
            if(tzm < tx0) answer|=1u;	// set bit at position 0
            return (int) answer;
        }
    }
    else
    {
        if(ty0 > tz0)
        { // PLANE XZ
            if(txm < ty0) answer|=4u;	// set bit at position 2
            if(tzm < ty0) answer|=1u;	// set bit at position 0
            return (int) answer;
        }
    }
    // PLANE XY
    if(txm < tz0) answer|=4u;	// set bit at position 2
    if(tym < tz0) answer|=2u;	// set bit at position 1
    return (int) answer;
}


int new_node(double txm, int x, double tym, int y, double tzm, int z)
{
    if(txm < tym)
    {
        if(txm < tzm){return x;}  // YZ plane
    }
    else
    {
        if(tym < tzm){return y;} // XZ plane
    }

    return z; // XY plane;
}


void proc_subtree(double tx0, double ty0, double tz0,
                  double tx1, double ty1, double tz1,
                  unsigned int depth, int justCreated,
                  Node* n, Node *parent, unsigned char childIndex, unsigned char a, Ray* r)
{
    double txm, tym, tzm;
    int currentNode;

    if (any_is_less(r->t_end, r->t_end, r->t_end, tx0, ty0, tz0))
    {
        // The ray endpoint happened before at least one of the minimum t values of this subtree, meaning the ray
        // ends before this subtree. Therefore, we don't want to do anything to this subtree.
        return;
    }

    if (tx1 < 0.0 || ty1 < 0.0 || tz1 < 0.0)
    {
        return;
    }
    else if (n == NULL && parent != NULL)
    {
        justCreated = createChildIfItDoesntExist(parent, childIndex);
        n = parent->children[childIndex];
    }

    if (!nodeHasAnyChildren(n))
    {
        if (depth == MAX_DEPTH)
        {
            double logLikelihoodUpdate = 0.0;

            if (any_is_greater(r->t_end, r->t_end, r->t_end, tx1, ty1, tz1))
            {
                // The ray endpoint does not occur in this subtree, but the ray passes through this subtree on its
                // way to the endpoint, and we're at our maximum depth. Therefore we need to give this node a vote
                // that it is free space.
                logLikelihoodUpdate = PROB_MISS_LOG;
            }
            else
            {
                // The ray endpoint occurs within this subtree, and we're at our maximum depth. Therefore we need to
                // give this node a vote that it is free space.
                logLikelihoodUpdate = PROB_HIT_LOG;
            }
            // Do the update
            n->logOdds += logLikelihoodUpdate;
            // Clamp the logOdds between the min/max
            n->logOdds = fmax(CLAMPING_THRES_MIN, fmin(n->logOdds, CLAMPING_THRES_MAX));
            return;
        }
        else if (!justCreated)
        {
            // Either the ray endpoint occurs within this subtree or the ray passes through this subtree on its way to
            // the endpoint, but we're not at our maximum depth yet and this is a previously pruned node, so we want to
            // expand this node and then proceed as normal.
            expandPrunedNode(n);
        }
    }

    txm = 0.5 * (tx0 + tx1);
    tym = 0.5 * (ty0 + ty1);
    tzm = 0.5 * (tz0 + tz1);

    currentNode = first_node(tx0, ty0, tz0, txm, tym, tzm);

    int createdChild = FALSE;
    do
    {
        switch (currentNode)
        {
            case 0:
                proc_subtree(tx0, ty0, tz0, txm, tym, tzm, depth + 1, createdChild, n->children[a], n, a, a, r);
                currentNode = new_node(txm, 4, tym, 2, tzm, 1);
                break;

            case 1:
                proc_subtree(tx0, ty0, tzm, txm, tym, tz1, depth + 1, createdChild, n->children[1u^a], n, 1u^a, a, r);
                currentNode = new_node(txm, 5, tym, 3, tz1, 8);
                break;

            case 2:
                proc_subtree(tx0, tym, tz0, txm, ty1, tzm, depth + 1, createdChild, n->children[2u^a], n, 2u^a, a, r);
                currentNode = new_node(txm, 6, ty1, 8, tzm, 3);
                break;

            case 3:
                proc_subtree(tx0, tym, tzm, txm, ty1, tz1, depth + 1, createdChild, n->children[3u^a], n, 3u^a, a, r);
                currentNode = new_node(txm, 7, ty1, 8, tz1, 8);
                break;

            case 4:
                proc_subtree(txm, ty0, tz0, tx1, tym, tzm, depth + 1, createdChild, n->children[4u^a], n, 4u^a, a, r);
                currentNode = new_node(tx1, 8, tym, 6, tzm, 5);
                break;

            case 5:
                proc_subtree(txm, ty0, tzm, tx1, tym, tz1, depth + 1, createdChild, n->children[5u^a], n, 5u^a, a, r);
                currentNode = new_node(tx1, 8, tym, 7, tz1, 8);
                break;

            case 6:
                proc_subtree(txm, tym, tz0, tx1, ty1, tzm, depth + 1, createdChild, n->children[6u^a], n, 6u^a, a, r);
                currentNode = new_node(tx1, 8, ty1, 8, tzm, 7);
                break;

            case 7:
                proc_subtree(txm, tym, tzm, tx1, ty1, tz1, depth + 1, createdChild, n->children[7u^a], n, 7u^a, a, r);
                currentNode = 8;
                break;

            default:
                assert(0);
        }
    } while (currentNode < 8);

    n->logOdds = maxChildLogLikelihood(n);
}


void ray_parameter(Octree* tree, Ray* r) {
    int createdRoot = FALSE;
    if (tree->root == NULL) {
        // Using calloc instead of malloc initializes the memory to zero, which means that the the root's `children`
        // array will be filled with zeros (which is equivalent to filling it will NULL pointers). It also sets the
        // root's log odds to 0, which we want as our initial value.
        tree->root = (Node *) calloc(1, sizeof(Node));
        createdRoot = TRUE;
    }

    unsigned char a = 0;

    // Since the tree is centered at (0, 0, 0), reflecting the origins is as simple as negating them.
    if (r->direction.x < 0.0f) {
        r->origin.x = -r->origin.x;
        r->direction.x = -r->direction.x;
        a |= 4u;
    }

    if (r->direction.y < 0.0f) {
        r->origin.y = -r->origin.y;
        r->direction.y = -r->direction.y;
        a |= 2u;
    }

    if (r->direction.z < 0.0f) {
        r->origin.z = -r->origin.z;
        r->direction.z = -r->direction.z;
        a |= 1u;
    }

    // Improve IEEE double stability
    double rdxInverse = 1.0 / r->direction.x;
    double rdyInverse = 1.0 / r->direction.y;
    double rdzInverse = 1.0 / r->direction.z;

    double tx0 = (tree->min.x - r->origin.x) * rdxInverse;
    double tx1 = (tree->max.x - r->origin.x) * rdxInverse;
    double ty0 = (tree->min.y - r->origin.y) * rdyInverse;
    double ty1 = (tree->max.y - r->origin.y) * rdyInverse;
    double tz0 = (tree->min.z - r->origin.z) * rdzInverse;
    double tz1 = (tree->max.z - r->origin.z) * rdzInverse;

    if (MAX(MAX(tx0, ty0), tz0) < MIN(MIN(tx1, ty1), tz1))
    {
        proc_subtree(tx0, ty0, tz0, tx1, ty1, tz1, 0, createdRoot, tree->root, NULL, 0, a, r);
    }
    else
    {
        printf("Ray outside of tree bounds\n");
    }
}

static inline void pruneNode(Node* node)
{
    if (nodeIsPrunable(node))
    {
        deleteAllChildren(node);
    }
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

    // for each point, create ray and call ray_parameter. After all points done, prune the tree
    for (size_t i = 0; i < numPoints; ++i)
    {
        Ray currentRay;
        initRay(&currentRay,
                sensorOrigin->x, sensorOrigin->y, sensorOrigin->z,
                points[i].x, points[i].y, points[i].z);

        ray_parameter(tree, &currentRay);
    }
}
