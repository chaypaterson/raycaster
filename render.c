#include <stdlib.h>
#include <string.h>

// This program renders images to PPM files.

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

void draw_frame(struct VoxelCube cube) {
    // Draw a white cage:
    for (unsigned row = 0; row < cube.resol.x; ++row) {
        for (char ch = 0; ch < 3; ++ch) {
            cube.buff[row][0][0][ch] = 100;
            cube.buff[row][cube.resol.y-1][0][ch] = 100;
            cube.buff[row][0][cube.resol.z-1][ch] = 100;
            cube.buff[row][cube.resol.y-1][cube.resol.z-1][ch] = 100;
        }
    }
    for (unsigned col = 0; col < cube.resol.y; ++col) {
        for (char ch = 0; ch < 3; ++ch) {
            cube.buff[0][col][0][ch] = 100;
            cube.buff[cube.resol.x-1][col][0][ch] = 100;
            cube.buff[0][col][cube.resol.z-1][ch] = 100;
            cube.buff[cube.resol.x-1][col][cube.resol.z-1][ch] = 100;
        }
    }
    for (unsigned lyr = 0; lyr < cube.resol.z; ++lyr) {
        for (char ch = 0; ch < 3; ++ch) {
            cube.buff[0][0][lyr][ch] = 100;
            cube.buff[cube.resol.x-1][0][lyr][ch] = 100;
            cube.buff[0][cube.resol.y-1][lyr][ch] = 100;
            cube.buff[cube.resol.x-1][cube.resol.y-1][lyr][ch] = 100;
        }
    }
}

void draw_rgb(struct VoxelCube cube) {
    // Draw an RGB cube in the voxel buffer:
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
}

void save_cube(struct VoxelCube cube, char *filename) {
    // Save a cube's voxel buffer to a file
    FILE *file = fopen(filename, "wb");
    const char signature[6] = "Voxel\n";
    fwrite(signature, sizeof(char), 6, file);

    // write resolution to file
    unsigned resol[3] = {cube.resol.x, cube.resol.y, cube.resol.z};
    fwrite(resol, sizeof(unsigned), 3, file);

    // write voxels to file
    for (unsigned row = 0; row < cube.resol.x; ++row) {
        for (unsigned col = 0; col < cube.resol.y; ++col) {
            for (unsigned lyr = 0; lyr < cube.resol.z; ++lyr) {
                fwrite(cube.buff[row][col][lyr], sizeof(float), 3, file);
            }
        }
    }

    fclose(file);
}

struct VoxelCube load_cube(char *filename) {
    // Read from a file into a voxel cube
    // The file format is: the first 6 characters are "Voxel\n", followed by
    // three unsigned integers: the dimensions of the cube. The rest of the file
    // consists of the contents of the voxel colour buffer.
    FILE *file = fopen(filename, "rb");
    // TODO 6 is a magic number here
    const char signature[6] = "Voxel\n";

    char header[6];
    fread(header, sizeof(char), 6, file);
    // check these two are equal:
    if (!strcmp(signature, header)) {
        printf("Wrong input file type for voxel cube\n");
        printf("%s vs %s\n", header, signature);
        exit(1);
    }

    unsigned resol[3]; // 3D resolution of cube
    fread(resol, sizeof(unsigned), 3, file);

    struct VoxelCube cube = new_unit_cube(resol[0], resol[1], resol[2]);

    // Load voxels:
    for (unsigned row = 0; row < cube.resol.x; ++row) {
        for (unsigned col = 0; col < cube.resol.y; ++col) {
            for (unsigned lyr = 0; lyr < cube.resol.z; ++lyr) {
                fread(cube.buff[row][col][lyr], sizeof(float), 3, file);
            }
        }
    }

    fclose(file);

    return cube;
}

int main() {
    // Construct voxel cube:
    struct VoxelCube cube = new_unit_cube(64, 64, 64);

    // Fill cube with stuff:
    draw_rgb(cube);
    draw_frame(cube);

    // Try saving and loading cube:
    printf("Saving cube...\n");
    save_cube(cube, "rgb.cube");
    free_unit_cube(cube); // flush the cube buffer
    printf("Loading cube...\n");
    cube = load_cube("rgb.cube"); // reload

    // unit test: check if the centre of the cube is inside the cube:
    printf("I should be 1: %d \n", is_inside_box(cube.geom.centre, cube));

    // construct imaging plane:
    struct ImagePlane plane = new_image_plane(640, 640);
    plane.geom.dims[0] = 2;
    plane.geom.dims[1] = 2;
    printf("Plane geometry: %g x %g\n", plane.geom.dims[0], plane.geom.dims[1]);

    // Set scene geometry:
    double theta = -M_PI * 0.25;
    double phi = 0.0f;

    // Create video with multiple views of the same cube:
    int maxframes = 150;

    for (int frame = 0; frame < maxframes; ++frame) {
        orient_image_plane(&plane, cube, 2.0, theta, phi);

        // render scene:
        raycast(plane, cube);

        // quantise colours and write out plane image buffer:
        char frameppm[sizeof("frame000.ppm")];
        sprintf(frameppm, "frame%03d.ppm", frame);

        FILE* img = fopen(frameppm, "w");

        putheader(img, plane.resol.cols, plane.resol.rows);

        for (unsigned row = 0; row < plane.resol.rows; ++row) {
            for (unsigned col = 0; col < plane.resol.cols; ++col) {
                Pixel colour;
                for (char ch = 0; ch < 3; ++ch) {
                    // Quantise with a saturation value and gamma correction:
                    colour[ch] = quantise(plane.buff[row][col][ch], 1.0, 0.75);
                    // wipe the image plane buffer:
                    plane.buff[row][col][ch] = 0.0f;
                }
                putpixel(img, colour);
            }
        }

        fclose(img);

        phi += 2.0 * M_PI / maxframes;
        theta += 0.25 * sin(frame * M_PI * 2.0 / maxframes) * 2.0 * M_PI/ maxframes;
    }

    free_unit_cube(cube);
    free_image_plane(plane);

    return 0;
}
