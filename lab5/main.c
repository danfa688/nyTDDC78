#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>
#include <math.h>
#include "physics.h"
#include "coordinate.h"
#include "extra.h"

int main (int argc, char *argv[]) {
	if(argc = 54){
		printf("usage: mpirun N<number_of_processors> program xprocesses yprocesses\n");
		exit(1);
	}
	int np, me;
	//Init MPI
	MPI_Init( &argc, &argv );
    MPI_Comm com = MPI_COMM_WORLD;
    MPI_Comm_size( com, &np );
    MPI_Comm_rank( com, &me );
	MPI_Status status;
	if ((atoi(argv[1]) * atoi(argv[2])) != np){
		printf("xprocesses * yprocess must be equal to number of processes allocated in mpirun!\n");
		MPI_Finalize();
		exit(1);
	}
	//Create mpi structure
    particle_t item;

	MPI_Datatype particle_mpi;
	int block_lengths [] = {1 , 1, 1, 1};
	MPI_Datatype block_types [] = {MPI_FLOAT, MPI_FLOAT, MPI_FLOAT, MPI_FLOAT, MPI_INT};
	MPI_Aint start, displ[5];

	MPI_Address(&item, &start);
	MPI_Address(&item.pcord.x, &displ[0]);
	MPI_Address(&item.pcord.y, &displ[1]);
	MPI_Address(&item.pcord.vx, &displ[2]);
	MPI_Address(&item.pcord.vy, &displ[3]);
	MPI_Address(&item.ptype, &displ[4]);

	displ[0] -= start;
	displ[1] -= start;
	displ[2] -= start;
	displ[3] -= start;
	displ[4] -= start;
	MPI_Type_struct(5, block_lengths, displ, block_types, &particle_mpi);

	MPI_Type_commit( &particle_mpi);
	//STOP Create mpi structure
	
	
	//Array with all particles
	particle_t* particles_array;
	particles_array = malloc(MAX_NO_PARTICLES*sizeof(*particles_array));
	
	//Randomize all particles position, velocity and starting angle
	int i, r, theta;
	for(i=0; i< INIT_NO_PARTICLES; i++)
	{
		r = rand()*MAX_INITIAL_VELOCITY;
		theta = rand()*2*PI;
		particles_array[i].pcord.vx = r*cos(theta);
		particles_array[i].pcord.vy = r*sin(theta);
		
		particles_array[i].pcord.x = rand()*BOX_HORIZ_SIZE;
		particles_array[i].pcord.y = rand()*BOX_VERT_SIZE;
	}
	
	// Loop for t seconds
	// Check particle collisions (collide)
	// Move particles part of collision (interact)
	// Move particles that has NOT been part of a collision (feuler)
	// Check wall collisions and add momentum (wall_collide)
	// End loop
	// Sum all momentum absorbed by the wall and divide this with WALL_LENGTH to get the pressure
	
	
}


