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

    // test access: fill buffer with zeroes
    for (unsigned ix = 0; ix < x; ++ix) {
        buff[ix] = malloc(sizeof(Colour[y][z]));
        for (unsigned iy = 0; iy < y; ++iy) {
            buff[ix][iy] = malloc(sizeof(Colour[z]));
            for (unsigned iz = 0; iz < z; ++iz) {
                buff[ix][iy][iz] = malloc(sizeof(Colour));
                for (char ch = 0; ch < 3; ++ch) {
                    buff[ix][iy][iz][ch] = 0;
                }
            }
        }
    }

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

    return;
}
