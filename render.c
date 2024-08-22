#include <stdlib.h>
#include <string.h>

// This program prints a binary PPM image to a file.

#include "pixel.h"
#include "raycast.h"

void gradient_test() {
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
}

int main() {
    // Test routines from raycast:
    struct VoxelCube cube = new_unit_cube(32, 32, 32);
    struct ImagePlane plane = new_image_plane(64, 64);

    // Set scene geometry:
    orient_image_plane(&plane, cube, 2.0, 1.0f, 1.0f);

    // Fill cube with stuff:
    for (unsigned row = 0; row < cube.resol.x; ++row) {
        for (unsigned col = 0; col < cube.resol.y; ++col) {
            for (unsigned lyr = 0; lyr < cube.resol.z; ++lyr) {
                for (char ch = 0; ch < 3; ++ch) {
                    // TODO: number of channels is usually 3 but this should
                    // really be a variable defined in a header file. what if we
                    // wanted an alpha channel?
                    cube.buff[row][col][lyr][ch] = 1.0f;
                }
            }
        }
    }

    // render:
    raycast(plane, cube);

    // quantise colours and write out plane image buffer:
    FILE* img = fopen("rendering1.ppm", "w");

    putheader(img, plane.resol.cols, plane.resol.rows);

    for (unsigned row = 0; row < plane.resol.rows; ++row) {
        for (unsigned col = 0; col < plane.resol.cols; ++col) {
            Pixel colour;
            for (char ch = 0; ch < 3; ++ch) {
                colour[ch] = quantise(plane.buff[row][col][ch]);
            }
            putpixel(img, colour);
        }
    }

    fclose(img);

    free_unit_cube(cube);
    free_image_plane(plane);

    return 0;
}
