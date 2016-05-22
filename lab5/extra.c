#include "extra.h"
#include <math.h>

//------------------------------ START OF INIT -------------------------------
void init_wall(cord_t* wall){
	wall->x0 = 0;
	wall->y0 = 0;
	wall->x1 = BOX_HORIZ_SIZE;
	wall->y1 = BOX_VERT_SIZE;
}

void create_my_area(area_t *my_a, int pid, int width, int height, int npx, int npy){
	// Calculating the x and y coordinates for the local area
	coord_calc(my_a, pid, width, height, npx, npy);
	// Calculating neighbours of the local area
	neigh_calc(my_a, pid, npx, npy);
}

void coord_calc(area_t* my_a, int pid, int width, int height, int npx, int npy){
	int x0, x1, y0, y1, local_size, tmp, rest;
	
	// Calculating X-coordinates
	local_size = width/npx;
	rest = width%npx;
	tmp = pid%npx;
	if(tmp < rest){
		x0 = tmp * (local_size + 1);
		x1 = x0 + local_size;
	}else if (tmp == rest){
		x0 = tmp * (local_size + 1);
		x1 = x0 + local_size - 1;
	}else{
		x0 = tmp * local_size + rest;
		x1 = x0 + local_size - 1;
	}
	
	// Calculating Y-coordinates
	local_size = height/npy;
	rest = height%npy;
	tmp = pid/npy;
	if(tmp < rest){
		y0 = tmp * (local_size + 1);
		y1 = y0 + local_size;
	}else if (tmp == rest){
		y0 = tmp * (local_size + 1);
		y1 = y0 + local_size - 1;
	}else{
		y0 = tmp * local_size + rest;
		y1 = y0 + local_size - 1;
	}
	
	my_a->x0 = x0;
	my_a->x1 = x1;
	my_a->y0 = y0;
	my_a->y1 = y1;
}

void neigh_calc(area_t* my_a, int pid, int npx, int npy){
	int neighbour[3][3];
	int i,k;
	// Fill neighbour list
	neighbour[0][0] = pid-npx-1;
	neighbour[0][1] = pid-npx;
	neighbour[0][2] = pid-npx+1;
	neighbour[1][0] = pid-1;
	neighbour[1][2] = pid+1;
	neighbour[2][0] = pid+npx-1;
	neighbour[2][1] = pid+npx;
	neighbour[2][2] = pid+npx+1;
	
	//Remove neighbouring walls 
	if(pid%npx==0){
		neighbour[0][0] = -1;
		neighbour[1][0] = -1;
		neighbour[2][0] = -1;
	}
	if(!((pid+1)%npx)){
		neighbour[0][2] = -1;
		neighbour[1][2] = -1;
		neighbour[2][2] = -1;
	}
	if((pid-npx)<0){
		neighbour[0][0] = -1;
		neighbour[0][1] = -1;
		neighbour[0][2] = -1;
	}
	if((pid+npx)>(npx*npy-1)){
		neighbour[2][0] = -1;
		neighbour[2][1] = -1;
		neighbour[2][2] = -1;
	}
	//Don't send to yourself!
	neighbour[1][1] = -1;
	
	//memcpy(&(my_a->neighbour_list[0]), &(neighbour[0]), 9*sizeof(neighbour[0]));
	for(i=0; i<3; i++){
		for(k=0; k<3; k++){
			(my_a->neighbour_list)[i][k].pid= neighbour[i][k];
			(my_a->neighbour_list)[i][k].send_buffer_length = 0;
		}
	}
}
// Initiatez all particles by randomizing all particles 
// position, velocity and starting angle:
void init_particles(area_t* local_area){
  int i;
  float r, theta;
	for(i=0; i< INIT_NO_PARTICLES; i++)
	{
		int local_width, local_height;
		local_width = (local_area->x1) - (local_area->x0) + 1;
		local_height = (local_area->y1) - (local_area->y0) + 1;
		r = (float)rand()/(float)RAND_MAX*(float)MAX_INITIAL_VELOCITY;
		theta = (float)rand()/(float)RAND_MAX*(float)2*PI;
		local_area->particle_array[i].pcord.vx = r*cos(theta);
		local_area->particle_array[i].pcord.vy = r*sin(theta);
		
		local_area->particle_array[i].pcord.x = (float)rand()/(float)RAND_MAX*(float)local_width;
		local_area->particle_array[i].pcord.y = (float)rand()/(float)RAND_MAX*(float)local_height;
	}
	local_area->no_particles = INIT_NO_PARTICLES;
	local_area->moment = 0;
}

void init_MPI(particle_t* item, int* me, int* np, MPI_Comm com, int* argc, char* argv[]){	
	//Init MPI
    MPI_Comm_size( com, np );
    MPI_Comm_rank( com, me );
	MPI_Status status;

	int block_lengths [] = {1 , 1, 1, 1};
	MPI_Datatype block_types [] = {MPI_FLOAT, MPI_FLOAT, MPI_FLOAT, MPI_FLOAT, MPI_INT};
	MPI_Aint start, displ[5];

	MPI_Address(item, &start);
	MPI_Address(&(item->pcord.x), &displ[0]);
	MPI_Address(&(item->pcord.y), &displ[1]);
	MPI_Address(&(item->pcord.vx), &displ[2]);
	MPI_Address(&(item->pcord.vy), &displ[3]);
	MPI_Address(&(item->ptype), &displ[4]);

	displ[0] -= start;
	displ[1] -= start;
	displ[2] -= start;
	displ[3] -= start;
	displ[4] -= start;
	MPI_Type_struct(5, block_lengths, displ, block_types, &particle_mpi);

	MPI_Type_commit( &particle_mpi);
}

//------------------------------ END OF INIT -------------------------------




void simulate(area_t* local_area, MPI_Comm com, cord_t* wall, int no_steps){
	int i;
	for(i=0; i<no_steps; i++){
		time_step(local_area, wall);
		communicate(local_area, com);
	}
}

void time_step(area_t* local_area, cord_t* wall){
  particle_t *p1, *p2;
	float t;
	int i, j, bool_collide;
	//for each particle
	for(i=0; i<local_area->no_particles-1; i++){
		bool_collide = 0;
		p1=&(local_area->particle_array[i]);
		//for each remaining particle
		for(j=i+1; j<local_area->no_particles; j++){
		  p2=&(local_area->particle_array[j]);
			//If two particles collided...
			t=collide(&(p1->pcord),&(p2->pcord));
			if(t != -1){
				//...update their position and velocities ...
				interact(&(local_area->particle_array[i].pcord),&(local_area->particle_array[j].pcord),t);
				//...and update total momentum of the two particles.				
				local_area->moment += wall_collide(&(p1->pcord),*wall);
				local_area->moment += wall_collide(&(p2->pcord),*wall);
				//We dont have to check the second particles position again so it's swapped with
				//the next element.
				swap_particle(&(local_area->particle_array[j]),&(local_area->particle_array[i+1]));
				bool_collide = 1;
				//Increments the iterator since we dont have to update that other particles position again.
				i++;
				break;
			}
		}
		if(!bool_collide){
			feuler(&(local_area->particle_array[i].pcord),STEP_SIZE);
			local_area->moment += wall_collide(&(p1->pcord),*wall);
		}
	}
}

void communicate(area_t* local_area, MPI_Comm com){
	add_to_send_buffer(local_area);
	int i,j;
	int flag = 1;
	for(i=0; i<3; i++){
		for(j=0; j<3; j++){
			if((local_area->neighbour_list[i][j].pid !=-1)){
				MPI_Isend(local_area->neighbour_list[i][j].send_buffer, 
							local_area->neighbour_list[i][j].send_buffer_length,
							particle_mpi, local_area->neighbour_list[i][j].pid, 
							local_area->neighbour_list[i][j].send_buffer_length, com, 
							&(local_area->neighbour_list[i][j].send_request));			
			}
		}
	}
	
	for(i=2; i>=0; i--){
	  for(j=2; j>=0; j--){
	    if((local_area->neighbour_list[i][j].pid !=-1)){
	      MPI_Probe(local_area->neighbour_list[i][j].pid, MPI_ANY_TAG, com,
			&(local_area->neighbour_list[i][j].status));
							
	      MPI_Recv(local_area->neighbour_list[i][j].receive_buffer, 
		       local_area->neighbour_list[i][j].status.MPI_TAG,
		       particle_mpi, local_area->neighbour_list[i][j].pid, 
		       MPI_ANY_TAG, com, MPI_STATUS_IGNORE);		
		       add_particles_from_buffer(local_area, &(local_area->neighbour_list[i][j]));
	    }
	  }
	}
	
	for(i=2; i>=0; i--){
	  for(j=2; j>=0; j--){
	    if((local_area->neighbour_list[i][j].pid !=-1)){
				MPI_Wait(&(local_area->neighbour_list[i][j].send_request), MPI_STATUS_IGNORE);
				local_area->neighbour_list[i][j].send_buffer_length = 0;
	    }
	  }
	}
}

float calculate_pressure(area_t* local_area, MPI_Comm com, int time){
	float global_moment;
	MPI_Reduce(&(local_area->moment), &global_moment, 1, MPI_FLOAT, MPI_SUM, 0, com);
	return ((global_moment/WALL_LENGTH)/time);
}

void swap_particle(particle_t* p1, particle_t* p2){
  	int tmp_ptype = p1->ptype;
    float tmp_x = p1->pcord.x;
    float tmp_y = p1->pcord.y;
    float tmp_vx = p1->pcord.vx;
    float tmp_vy = p1->pcord.vy;  
	
	p1->ptype = p2->ptype;
	p1->pcord.x = p2->pcord.x;
	p1->pcord.y = p2->pcord.y;
	p1->pcord.vx = p2->pcord.vx;
	p1->pcord.vy = p2->pcord.vy;

	
	p2->ptype = tmp_ptype;
	p2->pcord.x = tmp_x;
	p2->pcord.y = tmp_y;
	p2->pcord.vx = tmp_vx;
	p2->pcord.vy = tmp_vy;  
}

void add_to_send_buffer(area_t* local_area){
	int i;	
	for(i=0; i<(local_area->no_particles); i++){
		if(boundary_collide(&(local_area->particle_array[i]), local_area)){
			move(&(local_area->particle_array[i]),&(local_area->particle_array[local_area->no_particles-1]));
			(local_area->no_particles)--;
			i--;
		}
	}
}

int boundary_collide(particle_t* p, area_t* local_area){
	int rvalue; 
	int x0,x1,y0,y1;
	particle_t* ptr;
	x0 = local_area->x0;
	x1 = local_area->x1;
	y0 = local_area->y0;
	y1 = local_area->y1;
	rvalue = 1;
    if(p->pcord.x < x0){
		if(p->pcord.y < y0){
			//Left - Up
			ptr = &(local_area->neighbour_list[0][0].
					send_buffer[(local_area->neighbour_list[0][0].send_buffer_length)++]);
			move(ptr,p);
		}else if(p->pcord.y > y1){
			//Left - Down	
			ptr = &(local_area->neighbour_list[2][0].
					send_buffer[(local_area->neighbour_list[2][0].send_buffer_length)++]);
			move(ptr,p);	
		}else{
			//Left - Middle
			ptr = &(local_area->neighbour_list[1][0].
					send_buffer[(local_area->neighbour_list[1][0].send_buffer_length)++]);
			move(ptr,p);
		}	
    }else if(p->pcord.x > x1){
		if(p->pcord.y < y0){
			//Right - Up
			ptr = &(local_area->neighbour_list[0][2].
					send_buffer[(local_area->neighbour_list[0][2].send_buffer_length)++]);
			move(ptr,p);
		}else if(p->pcord.y > y1){
			//Right - Down
			ptr = &(local_area->neighbour_list[2][2].
					send_buffer[(local_area->neighbour_list[2][2].send_buffer_length)++]);
			move(ptr,p);		
		}else{
			//Right - Middle
			ptr = &(local_area->neighbour_list[1][2].
					send_buffer[(local_area->neighbour_list[1][2].send_buffer_length)++]);
			move(ptr,p);
		}
    }else{
		if(p->pcord.y < y0){
			//Middle - Up
			ptr = &(local_area->neighbour_list[0][1].
					send_buffer[(local_area->neighbour_list[0][1].send_buffer_length)++]);
			move(ptr,p);
		}else if(p->pcord.y > y1){
			//Middle - Down		
			ptr = &(local_area->neighbour_list[2][1].
					send_buffer[(local_area->neighbour_list[2][1].send_buffer_length)++]);
			move(ptr,p);
		}else{
			//Local area
			rvalue = 0;
		}
	}
	return rvalue;
}

// -------------- Utility -------------------
void move(particle_t* dest, particle_t* src){
	dest->pcord.x = src->pcord.x;
	dest->pcord.y = src->pcord.y;
	dest->pcord.vx = src->pcord.vx;
	dest->pcord.vy = src->pcord.vy;
	dest->ptype = src->ptype;
}

void add_particles_from_buffer(area_t* local_area, neighbour* from_neighbour){
	int i;
	for(i=0; i<(from_neighbour->status.MPI_TAG); i++){
		move(&(local_area->particle_array[local_area->no_particles]), &(from_neighbour->receive_buffer[i]));
		(local_area->no_particles)++;
	}
}








