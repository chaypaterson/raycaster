#include "raycast.h"

/* raycast.c
 * Implementations for the functions declared in raycast.h. This library
 * performs volume ray casting of a voxel cloud in a unit cube onto an image
 * plane.
 */

// Memory management: constructors and destructors for the bounding box and
// image plane

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
                for (char ch = 0; ch < 3; ++ch) {
                    buff[ix][iy][iz][ch] = 0;
                }
            }
        }
    }

    struct VoxelCube unitcube = {
        .resol.x = res_x,
        .resol.y = res_y,
        .resol.z = res_z,

        .geom.dims = {1, 1, 1},
        .geom.centre = {0.5, 0.5, 0.5},
        .geom.orient = {{1, 0, 0},
                        {0, 1, 0},
                        {0, 0, 1},
                        },

        .buff = buff
    };

    return unitcube;
}

void free_unit_cube(struct VoxelCube unitcube) {
    // Destroy the VoxelCube's buffer

    for (unsigned ix = 0; ix < unitcube.resol.x; ++ix) {
        for (unsigned iy = 0; iy < unitcube.resol.y; ++iy) {
            for (unsigned iz = 0; iz < unitcube.resol.z; ++iz) {
                free(unitcube.buff[ix][iy][iz]);
            }
            free(unitcube.buff[ix][iy]);
        }
        free(unitcube.buff[ix]);
    }

    free(unitcube.buff);
}

struct ImagePlane new_image_plane(unsigned rows, unsigned cols) {
    // Return an image plane with the specified resolution in a default location
    // and orientation in the scene.
    Colour **buff = malloc(sizeof(Colour[rows][cols]));

    // zero initialise image buffer:
    for (unsigned row = 0; row < rows; ++row) {
        buff[row] = malloc(sizeof(Colour[cols]));
        for (unsigned col = 0; col < cols; ++col) {
            buff[row][col] = malloc(Colour_size);
            for (char ch = 0; ch < 3; ++ch) {
                buff[row][col][ch] = 0;
            }
        }
    }

    struct ImagePlane image_plane = {
        .resol.rows = rows,
        .resol.cols = cols,

        // Default scene geometry:

        .geom.dims = {3, 3},
        .geom.centre = {0, 0, 0},
        .geom.tangent = {{0, 0, 0}, {0, 0, 0}},
        .geom.eye = {3, 0, 0},

        .buff = buff
    };

    return image_plane;
}

void free_image_plane(struct ImagePlane plane) {
    for (unsigned row = 0; row < plane.resol.rows; ++row) {
        for (unsigned col = 0; col < plane.resol.cols; ++col) {
            free(plane.buff[row][col]);
        }
        free(plane.buff[row]);
    }

    free(plane.buff);
}

// Geometry and ray casting: place imaging plane in the scene, detect
// collisions, and shoot rays.

void orient_image_plane(struct ImagePlane *image_plane,
                        struct VoxelCube ref_cube,
                        double range, double theta, double phi) {
    // Move the image plane so it is at a distance range from the centre of the reference
    // cube, and oriented so that the vector connecting the cube's centre with
    // the plane's centre is normal to the plane.

    // compute centre:
    Vector deltar = {range * sin(theta) * cos(phi),
                     range * sin(theta) * sin(phi),
                     range * cos(theta)};

    // Move image plane to new location:
    for (char axis = 0; axis < 3; ++axis) {
        image_plane->geom.centre[axis] = ref_cube.geom.centre[axis];
        image_plane->geom.centre[axis] += deltar[axis];
    }

    // Place eye at 4 * range:
    for (char axis = 0; axis < 3; ++axis) {
        image_plane->geom.eye[axis] = ref_cube.geom.centre[axis];
        image_plane->geom.eye[axis] += 4 * deltar[axis];
    }

    // compute tangent vectors:
    // The "x" tangent vector (points to the next row) is the theta unit vector
    // in polar coords

    image_plane->geom.tangent[0][0] = cos(theta) * cos(phi);
    image_plane->geom.tangent[0][1] = cos(theta) * sin(phi);
    image_plane->geom.tangent[0][2] = sin(theta);

    // The "y" tangent vector (points to the next col) is the phi unit vector
    // in polar coords

    image_plane->geom.tangent[1][0] = -sin(phi);
    image_plane->geom.tangent[1][1] = +cos(phi);
    image_plane->geom.tangent[1][2] = 0;

}

char is_inside_box(Vector point, struct VoxelCube box) {
    // check if point is inside box
    char test = 1; // Really this should be a bool: there is no advantage in
                   // declaring this as a bool, but it is a helpful hint to
                   // future engineers.

    for (char axis = 0; axis < 3; ++axis) {
        double diff = point[axis] - box.geom.centre[axis];
        test = test && (diff > -box.geom.dims[axis] * 0.5);
        test = test && (diff < +box.geom.dims[axis] * 0.5);
    }

    return test;
}

void shoot_ray(Colour restrict result,
               Vector start, Vector dir, struct VoxelCube cube, double tmax) {
    // shoot a ray from start to start+tmax*dir
    // accumulate the colours in the cube's buffer
    // write out the total colour to the pointer "result"
    // The result colour should be nonaliasing: we do not refer to this pointer
    // by any other name in this function

    // choose a step: plausible step is cube dimensions / resolution
    double dt = cube.geom.dims[0] * 1.0 / cube.resol.x;

    double extinction = 1.00; // TESTING
    double weight = 1.0;

    // compute this once and use it later:
    Vector corner; // of the cube
    for (char axis = 0; axis < 3; ++axis) {
        corner[axis] = cube.geom.centre[axis] - 0.5 * cube.geom.dims[axis];
    }

    Vector ray; // shoot this ray

    for (double t = 0; t < tmax; t += dt) {
        for (char axis = 0; axis < 3; ++axis) {
            ray[axis] = start[axis];
            ray[axis] += t * dir[axis];
        }

        if (is_inside_box(ray, cube)) {
            // get the coords of the nearest voxel:

            Vector difference;
            for (char axis = 0; axis < 3; ++axis) {
                difference[axis] = ray[axis] - corner[axis];
                // normalise to range 0-1:
                difference[axis] /= cube.geom.dims[axis];
            }

            // round position in cube to nearest voxel:
            unsigned row = (0.5 + cube.resol.x * difference[0]);
            unsigned col = (0.5 + cube.resol.y * difference[1]);
            unsigned lyr = (0.5 + cube.resol.z * difference[2]);

            // Note: I prefer multiplying by tests here to avoid branching

            // bounds checking: ensure that none of these exceed the resolution
            row += (cube.resol.x - 1 - row) * (row >= cube.resol.x);
            col += (cube.resol.y - 1 - col) * (col >= cube.resol.y);
            lyr += (cube.resol.z - 1 - lyr) * (lyr >= cube.resol.z);

            // and are all positive
            row += -row * (row < 0);
            col += -col * (col < 0);
            lyr += -lyr * (lyr < 0);

            // get the colour of this voxel and update result:

            for (char ch = 0; ch < 3; ++ch) {
                result[ch] += cube.buff[row][col][lyr][ch] * dt * weight;
            }

            weight *= extinction;
        }
    }
}

#include <stdio.h> // DEBUG
void raycast(struct ImagePlane image_plane, struct VoxelCube cube) {
    // Fill the image_plane's image buffer by casting rays onto the cube from
    // each pixel
    // Note: we shouldn't need to pass pointers to the objects in this function
    // because the buffers are already indirected.

    // guess a tmax: twice the distance from the centre of the cube to the
    // centre of the plane.
    Vector normal;
    double dist;

    for (char axis = 0; axis < 3; ++axis) {
        normal[axis] = image_plane.geom.centre[axis] - cube.geom.centre[axis];
        dist += normal[axis] * normal[axis];
    }

    dist = sqrt(dist);
    double tmax = 2 * dist;

    // We will also want to use normal later:
    for (char axis = 0; axis < 3; ++axis) {
        normal[axis] *= -1.0 / dist;
    }

    // get the corner of the image plane, it will be useful later

    Vector corner;
    for (char axis = 0; axis < 3; ++axis) {
        corner[axis] = image_plane.geom.centre[axis];

        // go backwards from the centre along the plane tangent vectors:
        for (char tang_ax = 0; tang_ax < 2; ++tang_ax) {
            double delta = image_plane.geom.tangent[tang_ax][axis];
            // across half the width or height of the plane:
            delta *= image_plane.geom.dims[tang_ax];
            delta *= 0.5;
            corner[axis] -= delta;
        }
    }

    for (unsigned row = 0; row < image_plane.resol.rows; ++row) {
        for (unsigned col = 0; col < image_plane.resol.cols; ++col) {
            // get scene coordinates of this pixel in the image plane
            Vector ray;
            
            // start at the plane corner and slide along the tangent vectors:
            for (char axis = 0; axis < 3; ++axis) {
                ray[axis] = corner[axis];
                // first the x/row/theta direction:
                double delta = image_plane.geom.tangent[0][axis];
                delta *= image_plane.geom.dims[0];
                delta *= row * 1.0 / image_plane.resol.rows;
                ray[axis] += delta;

                // then the y/col/phi direction:
                delta = image_plane.geom.tangent[1][axis];
                delta *= image_plane.geom.dims[1];
                delta *= col * 1.0 / image_plane.resol.cols;
                ray[axis] += delta;
            }

            // ray is now at the pixel coordinates in the scene. Cast this ray
            // in the direction from the eye to the pixel.
            // Previous version: direction "normal" (to the plane -- towards the cube)

            Vector dirn;
            double norm = 0;
            for (char axis = 0; axis < 3; ++axis) {
                dirn[axis] = ray[axis] - image_plane.geom.eye[axis];
                norm += dirn[axis] * dirn[axis];
            }
            norm = sqrt(norm);
            for (char axis = 0; axis < 3; ++axis) dirn[axis] /= norm;

            Colour pixel = image_plane.buff[row][col];

            shoot_ray(pixel, ray, normal, cube, tmax);
        }
    }
}
