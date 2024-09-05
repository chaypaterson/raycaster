#include <stdlib.h>
#include <string.h>

// This program renders images to PPM files.

#include "pixel.h"
#include "raycast.h"
#include <time.h>

void draw_rgb(struct VoxelCube cube) {
    // Draw an RGB cube in the voxel buffer:
    for (unsigned row = 0; row < cube.resol.x; ++row) {
        for (unsigned col = 0; col < cube.resol.y; ++col) {
            for (unsigned lyr = 0; lyr < cube.resol.z; ++lyr) {
                // Make an RGB cube:
                float voxel[VChannels] = {row * 1.0f / cube.resol.x,
                                  col * 1.0f / cube.resol.y,
                                  lyr * 1.0f / cube.resol.z};

                memcpy(cube.block[row][col][lyr], voxel, Colour_size);
            }
        }
    }
}

void draw_rgb_axes(struct VoxelCube cube) {
    // paint axes onto a cube
    // Identify the row, col, and lyr axes in red, green, and blue
    for (unsigned row = 0; row < cube.resol.x; ++row) {
        cube.block[row][0][0][0] = 20.0f;
    }
    for (unsigned col = 0; col < cube.resol.y; ++col) {
        cube.block[0][col][0][1] = 20.0f;
    }
    for (unsigned lyr = 0; lyr < cube.resol.z; ++lyr) {
        cube.block[0][0][lyr][2] = 20.0f;
    }
}

struct VoxelCube rgb_test(char* filename) {
    // Default unit test: create a pretty RGB cube 
    // Construct voxel cube:
    struct VoxelCube cube = new_unit_cube(64, 64, 64);

    // Fill a default cube with stuff:
    draw_rgb(cube);
    draw_rgb_axes(cube);

    // Try saving and loading cube:
    printf("Saving cube...\n");
    save_cube(cube, filename);
    free_unit_cube(cube); // flush the cube buffer
    printf("Loading cube...\n");
    cube = load_cube(filename); // reload

    return cube;
}

int main(int argc, char* argv[]) {
    // Step 1: get a cube
    struct VoxelCube (*cube_get)(char* filename);

    // Default behaviour: draw an RGB cube
    cube_get = rgb_test;
    char* filename = "rgb.cube";

    // Optional argv behaviour: load a cube from a file. 
    // e.g.
    //      ./renderer --load [file.cube]
    const char* loadme = "--load";
    for (char* *arg = argv; *(arg + 1) != NULL; ++arg) {
        if (!strcmp(loadme, *arg)) {
            // Set the cube getter and filename:
            cube_get = load_cube;
            filename = *(arg + 1);
            printf("Loading %s...\n", filename);
        }
    }

    // Get the cube:
    struct VoxelCube cube = cube_get(filename);
    // Add some axes:
    draw_rgb_axes(cube);

    // sanity check:
    printf("Max colour value: %f\n", maximum_colour_value(cube));

    srand(0); // seed RNG for antialiasing

    // unit test: check if the centre of the cube is inside the cube:
    printf("I should be 1: %d \n", is_inside_box(cube.geom.centre, cube));

    // Step 2: draw the cube and make a video
    // construct imaging plane:
    //struct Camera camera = new_camera(1080, 1920); // UHD
    //camera.geom.dims[1] *= 1920 * 1.0 / 1080; // rescale camera
    struct Camera camera = new_camera(640, 640);
    printf("Plane geometry: %g x %g\n", camera.geom.dims[0], camera.geom.dims[1]);

    // Set scene geometry:
    double theta = -M_PI * 0.25;
    double phi = 0.0f;

    float exposure = 1.0f; // lower is brighter
    float gamma = 0.95; // lower is more compressed

    // Create video with multiple views of the same cube:
    int maxframes = 150;

    time_t t_start = time(NULL);
    printf("Rendering...\n");

    for (int frame = 0; frame < maxframes; ++frame) {
        // Place the camera in the scene:
        orient_camera(&camera, cube, 2.0, theta, phi);

        // render scene (take a picture):
        raycast(camera, cube);

        // Choose a filename for this image:
        char frameppm[sizeof("frame000.ppm")];
        sprintf(frameppm, "frame%03d.ppm", frame);

        // quantise colours and write out camera image buffer to file:
        //save_film(frameppm, camera, exposure, gamma);

        // clean the image camera so it is ready for the next frame:
        blank_film(camera);

        // Move the camera around in an interesting way:
        phi += 2.0 * M_PI / maxframes;
        theta += 0.25 * sin(frame * M_PI * 2.0 / maxframes) * 2.0 * M_PI/ maxframes;
    }

    free_unit_cube(cube);
    destroy_film(camera);

    // Save reel to disk:
    time_t t_end = time(NULL);
    printf("Rendering complete, %d frames in %d s (%gfps)\n",
           maxframes, t_end - t_start, maxframes * 1.0 / (t_end - t_start));

    printf("\nDone\n");

    return 0;
}
