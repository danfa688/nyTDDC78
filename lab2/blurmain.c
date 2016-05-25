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
pixel* intermediate;
struct thread_data{
	long threadId; 
	int radius;
	char *infile;
	char *outfile;
};
struct thread_data thread_data_array[NUM_THREADS];

pthread_mutex_t mutex;
pthread_cond_t cond;
long reading_file_done;
long n_t_done_filter_part_1;
long n_t_done_filter_part_2;

int lproblem[NUM_THREADS][2];
double w[MAX_RAD];
int xsize, ysize;

/* Function declaration */
void calculate_local_problem_size(const int x,const int y,const int np, int lproblem[][2]);
int write_txt (const char* fname, const int radius, const double time, const int np);

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
	intermediate = malloc(MAX_PIXELS*sizeof(*intermediate));
	int colmax;
	
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
	calculate_local_problem_size(xsize, ysize, NUM_THREADS, lproblem);
	
	pthread_mutex_lock(&mutex);
	reading_file_done = 1;
    pthread_cond_broadcast(&cond); 
    pthread_mutex_unlock(&mutex);

    printf("Calling filter\n");

    clock_gettime(CLOCK_REALTIME, &stime);
	
    blurfilter_part_1(xsize, src, radius, w, intermediate, lproblem[tId][0], lproblem[tId][1]);
    
    /* Check all t done filtering part 1 */
    pthread_mutex_lock(&mutex);
    n_t_done_filter_part_1++;
    pthread_cond_broadcast(&cond);
    while(n_t_done_filter_part_1<NUM_THREADS){
    	pthread_cond_wait(&cond, &mutex);
    }
    pthread_mutex_unlock(&mutex);
    
    blurfilter_part_2(xsize, ysize, src, radius, w, intermediate, lproblem[tId][0], lproblem[tId][1]);
    
    /* Check all t done filtering part 2 */
    pthread_mutex_lock(&mutex);
    n_t_done_filter_part_2++;
    pthread_cond_broadcast(&cond);
    while(n_t_done_filter_part_2<NUM_THREADS){
    	pthread_cond_wait(&cond, &mutex);
    }
    pthread_mutex_unlock(&mutex);
    	
    clock_gettime(CLOCK_REALTIME, &etime);

    printf("Filtering took: %g secs\n", (etime.tv_sec  - stime.tv_sec) +
	   1e-9*(etime.tv_nsec  - stime.tv_nsec)) ;

    /* write result */
    //printf("Writing output file\n");
    
    if(write_txt (outfile, radius, (etime.tv_sec  - stime.tv_sec) + 1e-9*(etime.tv_nsec  - stime.tv_nsec), NUM_THREADS) != 0)
      exit(1);
      
   return NULL;
}

void *other_t(void *tParam) {
	struct thread_data *myData;
	long tId;
	int radius;

	myData=(struct thread_data *) tParam;
	tId = myData->threadId;
	radius = myData->radius;

    pthread_mutex_lock(&mutex);
    while (!reading_file_done)
    {
    	pthread_cond_wait(&cond, &mutex);
    }
    pthread_mutex_unlock(&mutex);
    
    blurfilter_part_1(xsize, src, radius, w, intermediate, lproblem[tId][0], lproblem[tId][1]);
    
    /* Check all t done filtering part 1 */
    pthread_mutex_lock(&mutex);
    n_t_done_filter_part_1++;
    pthread_cond_broadcast(&cond);
    while(n_t_done_filter_part_1<NUM_THREADS){
    	pthread_cond_wait(&cond, &mutex);
    }
    pthread_mutex_unlock(&mutex);
    
    blurfilter_part_2(xsize, ysize, src, radius, w, intermediate, lproblem[tId][0], lproblem[tId][1]);
    
    /* Check all t done filtering part 2 */
    pthread_mutex_lock(&mutex);
    n_t_done_filter_part_2++;
    pthread_cond_broadcast(&cond);
    while(n_t_done_filter_part_2<NUM_THREADS){
    	pthread_cond_wait(&cond, &mutex);
    }
    pthread_mutex_unlock(&mutex);   
    
    return NULL;
}

int main (int argc, char ** argv) {
	/* Init */
	pthread_mutex_init(&mutex,NULL);
	pthread_cond_init(&cond, NULL);
	reading_file_done=0;
	n_t_done_filter_part_1=0;
	n_t_done_filter_part_2=0;

	int radius;
	char* imagename;
	char* filename;
	int counter, imagecounter;


	
	for(imagecounter = 0; imagecounter<4; imagecounter++)
    {
    	switch(imagecounter)
    	{
    		case 0: 
    		imagename = "im1.ppm";
			filename = "im1blur.txt";
			break;
			case 1: 
    		imagename = "im2.ppm";
			filename = "im2blur.txt";
			break;
			case 2: 
    		imagename = "im3.ppm";
			filename = "im3blur.txt";
			break;
			case 3: 
    		imagename = "im4.ppm";
			filename = "im4blur.txt";
			break;
    	}
    	
    	thread_data_array[0].infile=imagename;
	    thread_data_array[0].outfile=filename;
    	
		for(counter = 0; counter < 10; counter ++)
		{
			radius= pow(2,counter);

			if((radius > MAX_RAD) || (radius < 1)) {
				fprintf(stderr, "Radius (%d) must be greater than zero and less then %d\n", radius, MAX_RAD);
				exit(1);
			}
	
			/*Create treads*/
			pthread_t threads[NUM_THREADS];
			int ret;
			long t;
			for(t=0;t<NUM_THREADS;t++) {
				thread_data_array[t].threadId = t;
				thread_data_array[t].radius = radius;
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
		}
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

int write_txt (const char* fname, const int radius, const double time, const int np) {

  FILE * fp;
  int errno = 0;

  if (fname == NULL) fname = "\0";
  fp = fopen (fname, "a");
  if (fp == NULL) {
    fprintf (stderr, "write_txt failed to open %s: %s\n", fname,strerror (errno));
    return 1;
  }
  
  fprintf(fp, "%d %d %f \n", np, radius, time);
  if (fclose (fp) == EOF) {
    perror ("Close failed");
    return 3;
  }
  return 0;
}

