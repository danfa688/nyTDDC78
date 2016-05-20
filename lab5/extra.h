#ifndef _extra_h
#define _extra_h

#include "definitions.h" 


void init_wall(cord_t* wall);
int in_area(area_t *a, particle_t *p);

void create_my_area(area_t *a_list, int pid, int width, int height, int npx, int npy);

void coord_calc(area_t* my_a, int pid, int width, int height, int npx, int npy);

void neigh_calc(area_t* my_a, int pid, int npx, int npy);





#endif
