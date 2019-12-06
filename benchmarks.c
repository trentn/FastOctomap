#include <stdio.h>
#include <stdlib.h>
#include "FastOctree.h"
#include "fast_code_utils.h"


#define NUMPOINTS 251377

int main(){
    FILE *fp;
    fp = fopen("pointcloud.txt", "r");
    
    Vector3d origin;
    Vector3d endpoints[NUMPOINTS];
    

    char tmp[4];
    double tmp2[3];

    fscanf(fp, "%s %lf %lf %lf %lf %lf %lf", tmp, &tmp2[0], &tmp2[1], &tmp2[2], &(origin.x), &(origin.y), &(origin.z));
    
    for(int i = 0; i < NUMPOINTS; i++){
        fscanf(fp, "%lf %lf %lf", &((endpoints[i]).x), &((endpoints[i]).y), &((endpoints[i]).z));
    }

    Octree tree;
    initOctree(&tree);
    
    //unsigned long long start = rdtsc();
    insertPointCloud(&tree, endpoints, NUMPOINTS, &origin);
    //unsigned long long end = rdtsc();

    //printf("%llu\n", (end-start));

    fclose(fp);
}