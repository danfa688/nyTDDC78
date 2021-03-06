/*
  File: thresfilter.h

  Declaration of pixel structure and thresfilter function.
    
 */
#ifndef _THRESFILTER_H_
#define _THRESFILTER_H_
/* NOTE: This structure must not be padded! */
typedef struct _pixel {
    unsigned char r,g,b;
} pixel;

void thresfilter_part_1(const int xsize, const int ysize, pixel* src, int* sum, int ylength, int yline);

void thresfilter_part_2(const int xsize, const int ysize, pixel* src, const int global_sum, int ylength, int yline);

#endif
