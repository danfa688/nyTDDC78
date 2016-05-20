#include "extra.h"
#include <math.h>

void init_wall(cord_t* wall){
	wall.x0 = 0;
	wall.y0 = 0;
	wall.x1 = BOX_HORIZ_SIZE;
	wall.y1 = BOX_VERT_SIZE;
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
			(my_a->neighbour_list)[i][k]= neighbour[i][k];
		}
	}
}
// Initiatez all particles by randomizing all particles 
// position, velocity and starting angle:
void init_particles(area_t* local_area){
		int i, r, theta;
	for(i=0; i< INIT_NO_PARTICLES; i++)
	{
		int local_width, local_height;
		local_width = (local_area->x1) - (local_area->x0) + 1;
		local_height = (local_area->y1) - (local_area->y0) + 1;
		r = rand()*MAX_INITIAL_VELOCITY;
		theta = rand()*2*PI;
		local_area->particles_array[i].pcord.vx = r*cos(theta);
		local_area->particles_array[i].pcord.vy = r*sin(theta);
		
		local_area->particles_array[i].pcord.x = rand()*local_width;
		local_area->particles_array[i].pcord.y = rand()*local_height;
	}
}
void timestep(){
	particle_collision();
	feuler();
	wall_collision();
	communicate();
}























