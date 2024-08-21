#include <stdlib.h>
#include <string.h>

// This program prints a binary PPM image to a file.

#include "pixel.h"
#include "raycast.h"

int main() {
    // Test routines from raycast:
    struct VoxelCube cube = new_unit_cube(16, 16, 16);
    struct ImagePlane plane = new_image_plane(cube, 1, 0, 0, 16, 16);
    // TODO: struct arguments for plane constructor? one for coords one for
    // resolution?
    free_unit_cube(cube);
    free_image_plane(plane);

    // Write out a test image:
    const char* imagefile = "testimg.ppm";
    Pixel BLACK = {0, 0, 0};
    // NB: COLOUR_MAX, 255, and 0xFF are synonymous on most systems.
    Pixel MAGENTA = {COLOUR_MAX, 0x00, COLOUR_MAX};
    Pixel CYAN = {0x00, COLOUR_MAX, COLOUR_MAX};
    Pixel RED = {COLOUR_MAX, 0x00, 0x00};

    int cols = 640;
    int rows = 320;

    // Open a file and write an image to it:

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
