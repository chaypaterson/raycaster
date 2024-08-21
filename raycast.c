#include "raycast.h"
#include <stdio.h>

/* raycast.c
 * Implementations for the functions declared in raycast.h. This library
 * performs volume ray casting of a voxel cloud in a unit cube onto an image
 * plane.
 */

struct VoxelCube new_unit_cube(unsigned x, unsigned y, unsigned z) {
    // Create colour buffer with correct dimensions:
    Colour ***buff = malloc(sizeof(Colour[x][y][z]));

    // fill buffer with zeroes
    for (unsigned ix = 0; ix < x; ++ix) {
        buff[ix] = malloc(sizeof(Colour[y][z]));
        for (unsigned iy = 0; iy < y; ++iy) {
            buff[ix][iy] = malloc(sizeof(Colour[z]));
            for (unsigned iz = 0; iz < z; ++iz) {
                buff[ix][iy][iz] = malloc(Colour_size);
                for (char ch = 0; ch < 3; ++ch) {
                    buff[ix][iy][iz][ch] = 0;
                }
            }
        }
    }

    struct VoxelCube unitcube = {
        .resol.x = x,
        .resol.y = y,
        .resol.z = z,

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

struct ImagePlane new_image_plane(struct VoxelCube ref_cube, 
                                  double range, double theta, double phi,
                                  unsigned rows, unsigned cols) {
    // Return an image plane at distance range from the centre of the reference
    // cube, and oriented so that the vector connecting the cube's centre with
    // the plane's centre is normal to the plane.
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

        .geom.dims = {2, 2},
        .geom.centre = {0, 0, 0}, // determined by ref_cube and args
        .geom.tangent = {{0, 0, 0}, {0, 0, 0}}, // see above

        .buff = buff
    };

    // compute centre:
    // compute tangent vectors:

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
