#ifndef _extra_h
#define _extra_h

#include "definitions.h" 

struct area {
    float x0 ;
    float x1 ;
    float y0 ;
    float y1 ;
    int neighbour_list[8]
} ;



typedef struct area area_t ;


int in_area(area_t *a, particle_t *p);

void create_my_area(area_t *a_list, int pid);







#endif
