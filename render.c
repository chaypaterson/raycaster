#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>

// This program prints a binary PPM image to a file.

#define COLOUR_MAX UCHAR_MAX

void putheader(FILE *img, int cols, int rows) {
    fprintf(img, "P6\n");
    fprintf(img, "%d %d\n", cols, rows);
    fprintf(img, "%d\n", COLOUR_MAX);
}

typedef unsigned char Pixel[3];

void putpixel(FILE *img, Pixel pixel){
    for (char ch = 0; ch < 3; ++ch)
        putc(pixel[ch], img);
}

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

int main() {
    const char* imagefile = "testimg.ppm";
    Pixel BLACK = {0, 0, 0};
    // NB: COLOUR_MAX, 255, and 0xFF are synonymous on most systems.
    Pixel MAGENTA = {COLOUR_MAX, 0x00, COLOUR_MAX};
    Pixel CYAN = {0x00, COLOUR_MAX, COLOUR_MAX};
    Pixel RED = {COLOUR_MAX, 0x00, 0x00};

    int cols = 640;
    int rows = 320;

    FILE* img = fopen(imagefile, "w");

    putheader(img, cols, rows);

    // Draw the image:
    for (int row = 0; row < rows; ++row) {
        for (int col = 0; col < cols; ++col) {
            Pixel colour;
            memcpy(colour, MAGENTA, sizeof(colour));

            // Choose a colour:
            double f = col * 1.0/cols;
            mixcolours_y_corr(colour, f, CYAN);

            putpixel(img, colour);
        }
    }

    fclose(img);
    return 0;
}
