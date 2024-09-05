#include "raycast.h"

/* raycast.c
 * Implementations for the functions declared in raycast.h. This library
 * performs volume ray casting of a voxel cloud in a unit cube onto an image
 * plane.
 */

// Memory management: constructors and destructors for the bounding box and
// camera

struct VoxelCube new_unit_cube(unsigned res_x, unsigned res_y, unsigned res_z) {
    // Create colour buffer with correct dimensions:
    Colour ***buff = malloc(sizeof(Colour[res_x][res_y][res_z]));

    // fill buffer with zeroes
    for (unsigned ix = 0; ix < res_x; ++ix) {
        buff[ix] = malloc(sizeof(Colour[res_y][res_z]));
        for (unsigned iy = 0; iy < res_y; ++iy) {
            buff[ix][iy] = malloc(sizeof(Colour[res_z]));
            for (unsigned iz = 0; iz < res_z; ++iz) {
                buff[ix][iy][iz] = malloc(Colour_size);
                for (char ch = 0; ch < VChannels; ++ch) {
                    buff[ix][iy][iz][ch] = 0;
                }
            }
        }
    }

    struct VoxelCube unitcube = {
        .resol.x = res_x,
        .resol.y = res_y,
        .resol.z = res_z,

        // Default: unit cube centred on 1/2,1/2,1/2 and oriented in the obvious
        // way:
        .geom.dims = {1, 1, 1},
        .geom.centre = {0.5, 0.5, 0.5},
        .geom.orient = {{1, 0, 0},
                        {0, 1, 0},
                        {0, 0, 1},
                        },

        // Default: cube should be transparent
        .geom.extinction = 1.0,

        // Initial value: TODO this is not portable for non-unit cubes
        .geom.circumradius = 0.866,

        // Allocate memory for voxel buffer
        .block = buff
    };

    return unitcube;
}

void free_unit_cube(struct VoxelCube unitcube) {
    // Destroy the VoxelCube's buffer

    for (unsigned ix = 0; ix < unitcube.resol.x; ++ix) {
        for (unsigned iy = 0; iy < unitcube.resol.y; ++iy) {
            for (unsigned iz = 0; iz < unitcube.resol.z; ++iz) {
                free(unitcube.block[ix][iy][iz]);
            }
            free(unitcube.block[ix][iy]);
        }
        free(unitcube.block[ix]);
    }

    free(unitcube.block);
}

// file management: load and save voxel cube data
#include <stdio.h>
#include <string.h>
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
                fwrite(cube.block[row][col][lyr], sizeof(float), VChannels, file);
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

    char header[sizeof(signature)];
    fread(header, sizeof(char), sizeof(signature), file);
    // check these two are equal:
    if (!strcmp(signature, header)) {
        printf("Wrong input file type for voxel cube\n");
        printf("%s vs %s\n", header, signature);
        exit(1);
    }

    unsigned resol[3]; // 3D resolution of cube
    fread(resol, sizeof(unsigned), 3, file);
    printf("Resolution: %dx%dx%d\n", resol[0], resol[1], resol[2]);

    struct VoxelCube cube = new_unit_cube(resol[0], resol[1], resol[2]);

    // Load voxels:
    for (unsigned row = 0; row < cube.resol.x; ++row) {
        for (unsigned col = 0; col < cube.resol.y; ++col) {
            for (unsigned lyr = 0; lyr < cube.resol.z; ++lyr) {
                fread(cube.block[row][col][lyr], sizeof(float), VChannels, file);
            }
        }
    }

    fclose(file);

    return cube;
}

float maximum_colour_value(struct VoxelCube cube) {
    // check the largest colour value
    // check maximum value
    float maxval = 0;

    for (unsigned row = 0; row < cube.resol.x; ++row) {
        for (unsigned col = 0; col < cube.resol.y; ++col) {
            for (unsigned lyr = 0; lyr < cube.resol.z; ++lyr) {
                for (unsigned ch = 0; ch < VChannels; ++ch) {
                    if (cube.block[row][col][lyr][ch] > maxval) {
                        maxval = cube.block[row][col][lyr][ch];
                    }
                }
            }
        }
    }

    return maxval;
}

struct Camera new_camera(unsigned rows, unsigned cols) {
    // Return an camera with the specified resolution in a default location
    // and orientation in the scene.

    // zero initialise image buffer:
    float Black[VChannels] = {0, 0, 0};
    Colour **buff = malloc(sizeof(Colour[rows][cols]));
    for (unsigned row = 0; row < rows; ++row) {
        buff[row] = malloc(sizeof(Colour[cols]));
        for (unsigned col = 0; col < cols; ++col) {
            buff[row][col] = malloc(Colour_size);
            memcpy(buff[row][col], Black, Colour_size);
        }
    }

    struct Camera camera = {
        .resol.rows = rows,
        .resol.cols = cols,

        // Default scene geometry:

        .geom.dims = {2, 2},
        .geom.centre = {0, 0, 0},
        .geom.tangent = {{0, 1, 0}, {0, 0, 1}},
        .geom.normal = {1, 0, 0},
        .geom.eye = {1, 0, 0},

        .film = buff
    };

    return camera;
}

#include "pixel.h"
void save_film(char* frameppm, struct Camera plane,
               float exposure, float gamma) {
    // Write the image plane to a file, setting overall brightness (exposure)
    // and correcting for gamma.
    FILE* img = fopen(frameppm, "w");

    putheader(img, plane.resol.cols, plane.resol.rows);

    for (unsigned row = 0; row < plane.resol.rows; ++row) {
        for (unsigned col = 0; col < plane.resol.cols; ++col) {
            Pixel colour;
            for (char ch = 0; ch < VChannels; ++ch) {
                // Quantise with a saturation value and gamma correction:
                colour[ch] = quantise(plane.film[row][col][ch], exposure, gamma);
            }
            putpixel(img, colour);
        }
    }

    fclose(img);
}

void blank_film(struct Camera plane) {
    // wipe image film, filling it with empty pixels:
    float Black[VChannels] = {0, 0, 0};
    for (unsigned row = 0; row < plane.resol.rows; ++row) {
        for (unsigned col = 0; col < plane.resol.cols; ++col) {
            memcpy(plane.film[row][col], Black, Colour_size);
        }
    }
}

void destroy_film(struct Camera plane) {
    for (unsigned row = 0; row < plane.resol.rows; ++row) {
        for (unsigned col = 0; col < plane.resol.cols; ++col) {
            free(plane.film[row][col]);
        }
        free(plane.film[row]);
    }
    free(plane.film);
}

// Geometry and ray casting: place imaging plane in the scene, detect
// collisions, and shoot rays.

void orient_camera(struct Camera *camera,
                        struct VoxelCube ref_cube,
                        double range, double theta, double phi) {
    // Move the camera so it is at a distance range from the centre of the 
    // reference cube, and oriented so that the vector connecting the cube's
    // centre with the plane's centre is normal to the plane.

    // compute new normal:
    Vector normal = {sin(theta) * cos(phi),
                     sin(theta) * sin(phi),
                     cos(theta)};

    for (char axis = 0; axis < 3; ++axis) {
        camera->geom.normal[axis] = normal[axis];
    }

    // Move centre of camera to new location:
    for (char axis = 0; axis < 3; ++axis) {
        camera->geom.centre[axis] = ref_cube.geom.centre[axis];
        camera->geom.centre[axis] += range * normal[axis];
    }

    // Place eye at double the range from the cube centre:
    for (char axis = 0; axis < 3; ++axis) {
        camera->geom.eye[axis] = camera->geom.centre[axis];
        camera->geom.eye[axis] += range * normal[axis];
    }

    // compute tangent vectors:
    // The "x" tangent vector (points to the next row) is the theta unit vector
    // in polar coords

    camera->geom.tangent[0][0] = cos(theta) * cos(phi);
    camera->geom.tangent[0][1] = cos(theta) * sin(phi);
    camera->geom.tangent[0][2] = -sin(theta);

    // The "y" tangent vector (points to the next col) is the phi unit vector
    // in polar coords

    camera->geom.tangent[1][0] = -sin(phi);
    camera->geom.tangent[1][1] = +cos(phi);
    camera->geom.tangent[1][2] = 0;
}

_Bool is_inside_box(Vector point, struct VoxelCube box) {
    // check if point is inside box
    _Bool test = 1;

    for (char axis = 0; axis < 3; ++axis) {
        double diff = point[axis] - box.geom.centre[axis];
        test = test && (diff > -box.geom.dims[axis] * 0.5);
        test = test && (diff < +box.geom.dims[axis] * 0.5);
    }

    return test;
}

unsigned roundcoord(unsigned side_resolution, double difference) {
    // round a difference coordinate (i.e. relative to the bounding cube) to an
    // integer coordinate for the nearest voxel
    unsigned rounded = (0.5 + side_resolution * difference);
    // Note: I prefer multiplying by tests here to avoid branching. This is a
    // library, the compiler might not optimise branches away.

    // bounds checking: ensure that the result does not exceed the resolution
    // rounded = min(side_resolution - 1, rounded)
    rounded += (side_resolution - 1 - rounded) * (rounded >= side_resolution);
    // ensure positivity:
    // rounded = max(rounded, 0)
    rounded += -rounded * (rounded < 0);

    return rounded;
}

// some vector geometry routines to raycast.c

double dot(Vector a, Vector b) {
    double product = 0;
    for (char i = 0; i < 3; ++i) product += a[i] * b[i];
    return product;
}

double randdbl(double epsilon) {
    // return a random variate with zero mean and variance epsilon^2 / 4:
    return epsilon * rand() / RAND_MAX - epsilon * 0.5;
}

void randomise(Vector dirn, double epsilon) {
    // add tiny random offsets to the direction as an antialiasing measure
    for (char axis = 0; axis < 3; ++axis) {
        dirn[axis] += randdbl(epsilon);
    }
}

void shoot_ray(Colour restrict result,
               Vector start, Vector dir, struct VoxelCube cube) {
    // shoot a ray from start to a reasonable end point
    // accumulate the colours in the cube's buffer
    // write out the total colour to the pointer "result"
    // The result colour should be nonaliasing: we do not refer to this pointer
    // by any other name in this function

    // choose a step: plausible step is cube dimensions / resolution
    double dt = cube.geom.dims[0] * 1.0 / cube.resol.x;

    double transmission = 1.0; // use this to model extinction/opacity of the voxels

    // compute this once and use it later:
    Vector corner; // of the cube
    for (char axis = 0; axis < 3; ++axis) {
        corner[axis] = cube.geom.centre[axis] - 0.5 * cube.geom.dims[axis];
    }

    // Intelligently work out a good end point and skip useless iterations at
    // the start. To do this, pre-compute intersections with a sphere that
    // bounds the bounding cube (the circumsphere):
    Vector difference; // between pixel and cube centres
    for (char axis = 0; axis < 3; ++axis) {
        difference[axis] = start[axis] - cube.geom.centre[axis];
    }
    double b = dot(dir, difference);
    double c = dot(difference, difference);
    c -= cube.geom.circumradius * cube.geom.circumradius;

    if ((b > 0) || (b * b <= c)) return; // skip this pixel
    // otherwise:
    double twidth = sqrt(b * b - c);

    Vector ray; // shoot this ray
    double tstart = -b - twidth;
    double tstop = -b + twidth;
    // antialiasing: randomise the ray start position a little bit:
    tstart += randdbl(dt); // averages to zero

    for (double t = tstart; t < tstop; t += dt) {
        // compute position of tip of the ray:
        // ray = start + t * dir
        for (char axis = 0; axis < 3; ++axis) {
            ray[axis] = start[axis];
            ray[axis] += t * dir[axis];
        }

        if (is_inside_box(ray, cube)) {
            // convert a position in space into the coords of the nearest
            // voxel:

            Vector difference;
            for (char axis = 0; axis < 3; ++axis) {
                difference[axis] = ray[axis] - corner[axis];
                // normalise to range 0-1:
                difference[axis] /= cube.geom.dims[axis];
            }

            // Nearest neighbour interpolation:
            // round position in cube to nearest voxel:
            unsigned row = roundcoord(cube.resol.x, difference[0]);
            unsigned col = roundcoord(cube.resol.y, difference[1]);
            unsigned lyr = roundcoord(cube.resol.z, difference[2]);

            // get the colour of this voxel and update result:
            for (char ch = 0; ch < VChannels; ++ch) {
                result[ch] += cube.block[row][col][lyr][ch] * dt * transmission;
            }

            /* TODO voxels could have alpha channels instead? */
            transmission *= cube.geom.extinction;
        }
    }
}

void start_position(Vector ray, struct Camera plane, 
                    unsigned row, unsigned col) {
    // get in-plane pixel coordinates relative to the plane's centre:
    double x = (plane.resol.rows - row) * 1.0 / plane.resol.rows - 0.5;
    double y = (plane.resol.cols - col) * 1.0 / plane.resol.cols - 0.5;

    x *= plane.geom.dims[0];
    y *= plane.geom.dims[1];

    // move x * tangent[0] + y * tangent[1] away from the centre:
    for (char axis = 0; axis < 3; ++axis) {
        ray[axis] = plane.geom.centre[axis];
    }
    for (char axis = 0; axis < 3; ++axis) {
        ray[axis] += x * plane.geom.tangent[0][axis];
    }
    for (char axis = 0; axis < 3; ++axis) {
        ray[axis] += y * plane.geom.tangent[1][axis];
    }
}

void get_eye_direction(Vector dirn, struct Camera plane, Vector ray) {
    // initially set dirn = ray - eye position
    double norm = 0;
    for (char axis = 0; axis < 3; ++axis) {
        dirn[axis] = ray[axis] - plane.geom.eye[axis];
        norm += dirn[axis] * dirn[axis];
    }
    norm = sqrt(norm);
    // Normalise direction vector:
    for (char axis = 0; axis < 3; ++axis) 
        dirn[axis] /= norm;
}

void raycast(struct Camera camera, struct VoxelCube cube) {
    // Fill the camera's image buffer by casting rays onto the cube from
    // each pixel
    // Note: we shouldn't need to pass pointers to the objects in this function
    // because the buffers are already indirected.

    for (unsigned row = 0; row < camera.resol.rows; ++row) {
        for (unsigned col = 0; col < camera.resol.cols; ++col) {
            // set scene coordinates of this pixel in the camera
            Vector ray;
            start_position(ray, camera, row, col);
            
            // ray is now at the pixel coordinates in the scene. Cast this ray
            // in the direction "normal" (to the plane -- towards the cube)
            Vector dirn;
            for (char axis = 0; axis < 3; ++axis) 
                dirn[axis] = -camera.geom.normal[axis];

            // Assign a pixel from the image buffer:
            Colour pixel = camera.film[row][col];

            // increment pixel value by the brightness along the ray:
            shoot_ray(pixel, ray, dirn, cube);
        }
    }
}
