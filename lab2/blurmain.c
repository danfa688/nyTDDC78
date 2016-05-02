#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <pthread.h> 
#include "ppmio.h"
#include "blurfilter.h"
#include "gaussw.h"
#include <math.h>
#define NUM_THREADS 4
#define MAX_RAD 3000
 
/* Global variables*/
pixel* src;
struct thread_data{
	long threadId; 
	int radius;
	char *infile;
	char *outfile;
};
struct thread_data thread_data_array[NUM_THREADS];

pthread_mutex_t mutex;
pthread_cond_t cond;
long reading_file_done=0;

int lproblem[NUM_THREADS][2];
int ldata[NUM_THREADS][2];
double w[MAX_RAD];

/* Function declaration */
void calculate_local_problem_size(const int x,const int y,const int np, int lproblem[][2]);
void calculate_local_allocation_size(int ldata[][2], int lproblem[][2],int const ysize, const int np, const int radius);


void *root_t(void *tParam) {
	struct thread_data *myData;
	long tId;
	int radius;
	char *infile;
	char *outfile;
	
	myData=(struct thread_data *) tParam;
	tId = myData->threadId;
	radius = myData->radius;
	infile = myData->infile;
	outfile = myData->outfile;
	
 	struct timespec stime, etime;
	src = malloc(MAX_PIXELS*sizeof(*src));
	int xsize, ysize, colmax;
	
	/* read file */
	if(read_ppm (infile, &xsize, &ysize, &colmax, (char *) src) != 0){
		exit(1);
	}

	if (colmax > 255) {
		fprintf(stderr, "Too large maximum color-component value\n");
		exit(1);
	}

	printf("Has read the image, generating coefficients\n");
	
	get_gauss_weights(radius, w);

	pthread_mutex_lock(&mutex);
	calculate_local_problem_size(xsize, ysize, NUM_THREADS, lproblem);
	calculate_local_allocation_size(ldata, lproblem, ysize , NUM_THREADS, radius);
	
	reading_file_done = 1;
    pthread_cond_signal(&cond);
    
    pthread_mutex_unlock(&mutex);

    //printf("Calling filter\n");

    //clock_gettime(CLOCK_REALTIME, &stime);

    //blurfilter(xsize, ysize, src, radius, w);

    //clock_gettime(CLOCK_REALTIME, &etime);

    //printf("Filtering took: %g secs\n", (etime.tv_sec  - stime.tv_sec) +
	  // 1e-9*(etime.tv_nsec  - stime.tv_nsec)) ;

    /* write result */
    //printf("Writing output file\n");
    
    //if(write_ppm (argv[3], xsize, ysize, (char *)src) != 0)
      //exit(1);
    
}


void *other_t(void *tParam) {
	struct thread_data *myData;
	long tId;
	int radius;

	myData=(struct thread_data *) tParam;
	tId = myData->threadId;
	radius = myData->radius;

	double w[MAX_RAD];
    get_gauss_weights(radius, w);
    
    pthread_mutex_lock(&mutex);
    while (!reading_file_done)
    {
    	pthread_cond_wait(&cond, &mutex);
    }
    pthread_mutex_unlock(&mutex);
    
}


int main (int argc, char ** argv) {
	/* Init */
	pthread_mutex_init(&mutex,NULL);
	pthread_cond_init(&cond, NULL);

	int radius;
	/* Take care of arguments */
	if (argc != 4) {
		fprintf(stderr, "Usage: %s radius infile outfile\n", argv[0]);
		exit(1);
	}
	radius=atoi(argv[1]);
		if((radius > MAX_RAD) || (radius < 1)) {
		fprintf(stderr, "Radius (%d) must be greater than zero and less then %d\n", radius, MAX_RAD);
		exit(1);
	}
	thread_data_array[0].infile=argv[2];
	thread_data_array[0].outfile=argv[3];
	
	/*Create treads*/
	pthread_t threads[NUM_THREADS];
	int ret;
	long t;
	for(t=0;t<NUM_THREADS;t++) {
		thread_data_array[t].threadId = t;
		thread_data_array[t].radius = radius;
		printf("In main: creating thread %ld\n", t);
		if (t == 0) {
			ret = pthread_create(&threads[t], NULL, root_t, (void *)&thread_data_array[t]);
		}else{
			ret = pthread_create(&threads[t], NULL, other_t, (void *)&thread_data_array[t]);
		}
		if (ret) {
			printf("ERROR! Return code: %d\n", ret);
			exit(-1);
		}
	} 
	for(t=0;t<NUM_THREADS;t++) {
		pthread_join(threads[t], NULL);
	}
	
	return(0);
}



void calculate_local_problem_size(const int x,const int y,const int np, int lproblem[][2]){
	int linesize, rest;
	int lysize[np];
	linesize = floor(y/np);
	rest = y%np;
	int i;
	int k;
	for (i=0; i<np; i++){
		if(i<rest){
			lysize[i]=linesize+1;
		}else{
			lysize[i]=linesize;
		}
	}
	int tmp;
	for (i=0; i<np; i++){
		tmp=0;
		for (k=0; k<i; k++){
			tmp += lysize[k];
		}
		
		lproblem[i][0]=tmp;
		lproblem[i][1]=lysize[i];
	}
}

void calculate_local_allocation_size(int ldata[][2], int lproblem[][2],int const ysize , const int np, const int radius){
	int i;
	for(i=0; i<np; i++){
		if(lproblem[i][0] > radius){
			ldata[i][0]=lproblem[i][0]-radius;
		}else{
			ldata[i][0]=0;
		}
		
		if(lproblem[i][0]+lproblem[i][1]-1 + radius < ysize){
			ldata[i][1]=lproblem[i][0]-ldata[i][0]+lproblem[i][1] + radius;
		}else{
			ldata[i][1]=ysize-ldata[i][0];
		}
	}
}

