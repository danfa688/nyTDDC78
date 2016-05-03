/*
  File: blurfilter.h

  Declaration of pixel structure and blurfilter function.
    
 */

#ifndef _BLURFILTER_H_
#define _BLURFILTER_H_

/* NOTE: This structure must not be padded! */
typedef struct _pixel {
    unsigned char r,g,b;
} pixel;

void blurfilter_part_1(const int xsize, pixel* src, const int radius, const double *w, pixel* dst, const int startline, const int ylength);

void blurfilter_part_2(const int xsize, const int ysize, pixel* src, const int radius, const double *w, pixel* dst, const int startline, const int ylength);

#endif
