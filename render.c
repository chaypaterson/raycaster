#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// This program prints a binary PPM image to a file.

static int COLOUR_MAX = 255;

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

int main() {
    const char* imagefile = "testimg.ppm";
    Pixel BLACK = {0, 0, 0};
    Pixel MAGENTA = {255, 0x00, 255};
    Pixel CYAN = {0x00, 255, 255};
    Pixel RED = {255, 0x00, 0x00};

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
            mixcolours(colour, f, CYAN);

            putpixel(img, colour);
        }
    }

    fclose(img);
    return 0;
}
