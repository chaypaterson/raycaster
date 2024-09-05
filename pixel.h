#ifndef PIXEL_H
#define PIXEL_H

#include <stdio.h>
#include <limits.h>

/* Assume 8-bit colours, represented as unsigned characters.
 * Define a maximum value in a system-independent way.
 * A type to store 24-bit RGB pixels as three channels.
 */

#define COLOUR_MAX UCHAR_MAX
//const static unsigned COLOUR_MAX = UCHAR_MAX;
typedef unsigned char Pixel[3];

// File operations:
// Write a valid PPM header to a file:
void putheader(FILE *img, int cols, int rows);

// Write a single pixel to a file:
void putpixel(FILE *img, Pixel pixel);

// Pixel operations:
// Blend two colours:
void mixcolours(Pixel newcolour, double weight1, Pixel colour2);
// Compute luminance:
double luma_y(Pixel colour);
// Blend two colours, correcting for constant luminance:
void mixcolours_y_corr(Pixel newcolour, double weight1, Pixel colour2);

// Floating point colour conversion functions and helpers:

// quantise a floating point value to a char after applying gamma correction:
char quantise(float colourch, float sat, float gamma);

#endif //PIXEL_H
