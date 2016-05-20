#include "physics.h"
#include "coordinate.h"
#include "extra.h"

int main (int argc, char *argv[]) {

	int np, me, npx, npy;
	MPI_Comm com;
	MPI_Status status;
	particle_t item;
	init_MPI(&item, &np, &me, &com ,&argc, &argv);

	if(argc != 3){
		if(me == 0){
			fprintf(stderr,"usage: mpirun n<number_of_processors> program xprocesses yprocesses\n");
		}
		MPI_Finalize();			
		exit(1);
	}
	npx = atoi(argv[1]);
	npy = atoi(argv[2]); 
	
	if ((npx * npy) != np){
		if(me == 0){
			printf("xprocesses * yprocess must be equal to number of processes allocated in mpirun!\n");
		}
		MPI_Finalize();
		exit(1);
	}
	
	//Initiate wall - Same for every process
	cord_t wall;
	init_wall(&wall);
	//Create local area which each process is responsible over
	area_t local_area;
	create_my_area(&local_area, me, BOX_HORIZ_SIZE, BOX_VERT_SIZE, npx, npy);

	//Allocates space to all local particles and initiates them
	init_particles(&local_area);

	time_step(&local_area, &wall);
	// Loop for t seconds


	
	// Sum all momentum absorbed by the wall and divide this with WALL_LENGTH to get the pressure
	
	MPI_Finalize();
}


