#include "physics.h"
#include "coordinate.h"
#include "extra.h"

int main (int argc, char *argv[]) {

	int np, me, npx, npy, time;
	float pressure;
	MPI_Status status;
	particle_t item;
	MPI_Init( &argc, &argv );
	MPI_Comm com = MPI_COMM_WORLD;
	//	init_MPI(&item, &np, &me, com ,&argc, &argv);
	//Init MPI
    MPI_Comm_size( com, &np );
    MPI_Comm_rank( com, &me );

	int block_lengths [] = {1 , 1, 1, 1};
	MPI_Datatype block_types [] = {MPI_FLOAT, MPI_FLOAT, MPI_FLOAT, MPI_FLOAT, MPI_INT};
	MPI_Aint start, displ[5];

	MPI_Address(&item, &start);
	MPI_Address(&(item.pcord.x), &displ[0]);
	MPI_Address(&(item.pcord.y), &displ[1]);
	MPI_Address(&(item.pcord.vx), &displ[2]);
	MPI_Address(&(item.pcord.vy), &displ[3]);
	MPI_Address(&(item.ptype), &displ[4]);

	displ[0] -= start;
	displ[1] -= start;
	displ[2] -= start;
	displ[3] -= start;
	displ[4] -= start;
	MPI_Type_struct(5, block_lengths, displ, block_types, &particle_mpi);

	MPI_Type_commit( &particle_mpi);

	if(argc != 4){
		if(me == 0){
			fprintf(stderr,"usage: mpirun n<number_of_processors> program xprocesses yprocesses time\n");
		}
		MPI_Finalize();			
		exit(1);
	}
	npx = atoi(argv[1]);
	npy = atoi(argv[2]); 
	time = atoi(argv[3]);
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

	simulate(&local_area, com, &wall, time);
	pressure = calculate_pressure(&local_area, com, time);
	if(me==0){
		printf("Simulated pressure: %f \n", pressure);
		int i;
	printf("Pressure according to the ideal gas law: %f \n", (float)(2*INIT_NO_PARTICLES*np) /(float)(3*BOX_HORIZ_SIZE*BOX_VERT_SIZE)*0.5*1*powf(35,2));
	}
	
	MPI_Finalize();
}


