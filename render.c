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
    // Construct voxel cube:
    struct VoxelCube cube = new_unit_cube(32, 32, 32);

    // Fill cube with stuff:
    for (unsigned row = 0; row < cube.resol.x; ++row) {
        for (unsigned col = 0; col < cube.resol.y; ++col) {
            for (unsigned lyr = 0; lyr < cube.resol.z; ++lyr) {
                // Make an RGB cube:
                float voxel[3] = {row * 1.0f / cube.resol.x,
                                  col * 1.0f / cube.resol.y,
                                  lyr * 1.0f / cube.resol.z};

                for (char ch = 0; ch < 3; ++ch) {
                    // TODO: number of channels is usually 3 but this should
                    // really be a variable defined in a header file. what if we
                    // wanted an alpha channel?
                    cube.buff[row][col][lyr][ch] = voxel[ch];
                }
            }
        }
    }

    // unit test: check if the centre of the cube is inside the cube:
    printf("I should be 1: %d \n", is_inside_box(cube.geom.centre, cube));

    // construct imaging plane:
    struct ImagePlane plane = new_image_plane(640, 640);
    //plane.geom.dims[0] = 2;
    //plane.geom.dims[1] = 2;
    printf("Plane geometry: %g x %g\n", plane.geom.dims[0], plane.geom.dims[1]);

    // Set scene geometry:
    double theta = M_PI * 0.25;//39 * M_PI / 150;//M_PI * 0.25;
    double phi = 0.0f;

    // Create video with multiple views of the same cube:
    int maxframes = 150;

    // TODO DEBUG: perspective looks weird at certain angles theta, clearly
    // there's a mistake with geometry somewhere

    for (int frame = 0; frame < maxframes; ++frame) {
        orient_image_plane(&plane, cube, 2.0, theta, phi);

        // render scene:
        raycast(plane, cube);

        // quantise colours and write out plane image buffer:
        char frameppm[sizeof("frame000.ppm")];
        sprintf(frameppm, "frame%03d.ppm", frame);
        //FILE* img = fopen("rendering1.ppm", "w");
        FILE* img = fopen(frameppm, "w");

        putheader(img, plane.resol.cols, plane.resol.rows);

        for (unsigned row = 0; row < plane.resol.rows; ++row) {
            for (unsigned col = 0; col < plane.resol.cols; ++col) {
                // DEBUG:
                if (row == 0 && col == 0.5 * plane.resol.cols) {
                    printf("intended normal:\n");
                    Vector normal = {-sin(theta) * cos(phi), -sin(theta) * sin(phi), -cos(theta)};
                    printf("[%g %g %g]\n", normal[0], normal[1], normal[2]);
                }

                Pixel colour;
                for (char ch = 0; ch < 3; ++ch) {
                    colour[ch] = quantise(plane.buff[row][col][ch], 1.0, 0.75);
                    // wipe the image plane buffer:
                    plane.buff[row][col][ch] = 0.0f;
                }
                putpixel(img, colour);
            }
        }

        fclose(img);

        //phi += 2.0 * M_PI / maxframes;
        theta -= M_PI * 1.0 / maxframes;
    }

    free_unit_cube(cube);
    free_image_plane(plane);

    return 0;
}
