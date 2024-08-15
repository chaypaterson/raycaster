#include <stdio.h>
#include <stdlib.h>

// This program prints a binary PPM image to a file.

static int COLOUR_MAX = 255;

void putheader(FILE *img, int cols, int rows) {
    fprintf(img, "P6\n");
    fprintf(img, "%d %d\n", cols, rows);
    fprintf(img, "%d\n", COLOUR_MAX);
}

typedef char Pixel[3];

void putpixel(FILE *img, Pixel pixel){
    for (char ch = 0; ch < 3; ++ch)
        putc(pixel[ch], img);
}

int main() {
    const char* imagefile = "testimg.ppm";
    Pixel BLACK = {0, 0, 0};
    int cols = 640;
    int rows = 320;

    FILE* img = fopen(imagefile, "w");

    putheader(img, cols, rows);

    for (int row = 0; row < rows; ++row)
        for (int col = 0; col < cols; ++col)
            putpixel(img, BLACK);

    fclose(img);
    return 0;
}
