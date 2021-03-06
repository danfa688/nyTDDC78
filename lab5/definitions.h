#include<stdlib.h>
#include<math.h>

#include "coordinate.h"
#include "physics.h"

#ifndef _definitions_h
#define _definitions_h

#define PI 3.141592653

#define MAX_NO_PARTICLES  100000  /* Maximum number of particles/processor */
#define INIT_NO_PARTICLES 10000    /* Initial number of particles/processor */
#define MAX_INITIAL_VELOCITY 50


#define BOX_HORIZ_SIZE 10000.0
#define BOX_VERT_SIZE 10000.0
#define WALL_LENGTH (2.0*BOX_HORIZ_SIZE+2.0*BOX_VERT_SIZE)

#define PARTICLE_BUFFER_SIZE MAX_NO_PARTICLES/5
#define COMM_BUFFER_SIZE  5*PARTICLE_BUFFER_SIZE

typedef struct {
  pcord_t  pcord;
  int ptype;        /* Used to simulate mixing of gases */ 
} particle_t;

typedef struct{
	int pid;
	int send_buffer_length;
	int receive_buffer_length;
	particle_t* send_buffer;
	particle_t* receive_buffer;
	MPI_Status status;
	MPI_Request send_request;
	MPI_Request recv_request;
} neighbour;

typedef struct {
    float x0;
    float x1;
    float y0;
    float y1;
    neighbour neighbour_list[3][3];

	particle_t* particle_array;
	int no_particles;
	float moment;
} area_t;



#endif
