#include <stdlib.h>
#include <stdio.h>

typedef struct {
    float x0 ;
    float x1 ;
    float y0 ;
    float y1 ;
    int neighbour_list[3][3]
} area_t;

void neigh_calc(area_t* my_a, int pid, int npx, int npy){
	int neighbour[3][3];
	//neighbour = malloc(sizeof(int)*9);
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
	memcpy(&(my_a->neighbour_list[0]), &(neighbour[0]), 9*sizeof(neighbour[0]));
	/*for(i = 0; i < 3; i++){
		memcpy(&(my_a->neighbour_list[i]), &(neighbour[i*3]), 3*sizeof(neighbour[0]));
	}
	/*for(i=0; i<3; i++){
		for(k=0; k<3; k++){
			(my_a->neighbour_list)[i][k]= neighbour[i][k];
		}
	}*/
}

int main (void){
	int i,j;
	int pid = 4;
	int npx = 3;
	int npy = 3;
	
	area_t testArray;
    
	neigh_calc(&testArray,pid,npx,npy);
	
	/*for (i=0; i<3; i++){
		for(j=0; j<3; j++){*/
			printf("[%d]",testArray.neighbour_list[0][0]);
	/*	}
		printf("\n");
	}*/
	return 0;
}
