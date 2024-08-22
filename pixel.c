#include "pixel.h"

/* Library containing helper functions for PPM files and pixel operations.
 */

// File operations:

void putheader(FILE *img, int cols, int rows) {
    fprintf(img, "P6\n");
    fprintf(img, "%d %d\n", cols, rows);
    fprintf(img, "%d\n", COLOUR_MAX);
}

void putpixel(FILE *img, Pixel pixel) {
    for (char ch = 0; ch < 3; ++ch)
        putc(pixel[ch], img);
}

// Pixel operations:

void mixcolours(Pixel newcolour, double weight1, Pixel colour2) {
    // mix colour2 into newcolour
    // weight1 should be 0 < weight1 < 1
    for (char ch = 0; ch < 3; ++ch) {
        newcolour[ch] += (char)(weight1 * (colour2[ch] - newcolour[ch]));
    }
}

double luma_y(Pixel colour) {
    // return the sRGB luminance of a pixel
    //return 0.2126 * colour[0] + 0.7152 * colour[1] + 0.0722 * colour[2];
    return 1.0f * colour[0] + 1.0f * colour[1] + 1.0f * colour[2];
    //return colour[0] * colour[0] + colour[1] * colour[1] + colour[2] * colour[2];
}

void mixcolours_y_corr(Pixel newcolour, double weight1, Pixel colour2) {
    // mix the colours, but correct so that the luminance is constant before and
    // after:
    double y_before = luma_y(newcolour);
    mixcolours(newcolour, weight1, colour2);
    double y_after = luma_y(newcolour);

    for (char ch = 0; ch < 3; ++ch) {
        newcolour[ch] += (char)((y_before/y_after - 1.0) * newcolour[ch]);
    }
}

char quantise(float colourch) {
    // quantise a floating point colour channel value
    return (char)(colourch);
}
