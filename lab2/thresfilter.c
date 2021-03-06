#include "thresfilter.h"
#include <stdio.h>

void thresfilter_part_1(const int xsize, const int ysize, pixel* src, int* sum, int ylength, int yline){
#define uint unsigned int 

  uint lsum, i, psum, nump;

  nump = xsize * ylength;

  for( i = yline*xsize, lsum = 0; i < (nump + yline*xsize); i++) {
    lsum += (uint)src[i].r + (uint)src[i].g + (uint)src[i].b;
  }

  lsum /= nump;
  *sum = lsum;
}

void thresfilter_part_2(const int xsize, const int ysize, pixel* src, const int global_sum, int ylength, int yline){
#define uint unsigned int 

  uint i, psum, nump;

  nump = xsize * ylength;
  for(i = yline*xsize; i < (nump + yline*xsize); i++) {
    psum = (uint)src[i].r + (uint)src[i].g + (uint)src[i].b;
    if(global_sum > psum) {
      src[i].r = src[i].g = src[i].b = 0;
    }
    else {
      src[i].r = src[i].g = src[i].b = 255;
    }
  }
}


