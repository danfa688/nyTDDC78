#ifndef _extra_h
#define _extra_h

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <mpi.h>

#include "definitions.h" 



void init_wall(cord_t* wall);

int in_area(area_t *a, particle_t *p);

void create_my_area(area_t *a_list, int pid, int width, int height, int npx, int npy);

void coord_calc(area_t* my_a, int pid, int width, int height, int npx, int npy);

void neigh_calc(area_t* my_a, int pid, int npx, int npy);

void Init_MPI(particle_t* item,int* me, int* np, MPI_Comm* com, int* argc, char** argv);

void swap_particle(particle_t* p1, particle_t* p2);

void time_step(area_t* local_area, cord_t* wall);

int boundary_collide(particle_t* p, area_t* local_area);

void move(particle_t* dest, particle_t* src);
#endif
