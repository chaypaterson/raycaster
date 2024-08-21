#ifndef RAYCAST_H
#define RAYCAST_H
/* raycast.h
 * This library generates an image plane of a unit cube filled with voxels.
 * The image plane is a 2d array of colours: the internal representation of
 * colours in this library is as a tuple of 3 floats. These need to be quantised
 * down to 24-bit RGB before being written out to an actual image, because most
 * existing display hardware does not support true high dynamic range.
 *
 * The voxel cube is a 3D array of colours. It has an associated resolution
 * (width, length, and breadth), and physical dimensions (1 unit, 1 unit, 1
 * unit). The unit cube is centred on [1/2, 1/2, 1/2], not [0, 0, 0], the
 * origin. This is to simplify bounding box checks.
 *
 * Similarly, the image plane is a 2d array of colours. It has an associated
 * resolution (rows, columns) and physical dimensions (2 units by 2 units: any
 * view of the cube should fit inside with a little bit to spare).
 *
 * As well as the data stored in image planes and voxels, the cube and image
 * plane have geometry. They have location in a "physical" space. Points in
 * space are represented with vectors, which are arrays of doubles. We may as
 * well use double precision for the vectors, because we will never have to
 * remember very many at any given time, and will mainly allocate them on the
 * stack.
 *
 * The geometry of the cube and image plane are only completely specified with
 * orientations. The cube is assumed to be oriented along the XYZ axes in the
 * obvious way, and the orientation of the image plane is specified relative to
 * the cube. The orientation of the image plane is given with two tangent
 * vectors. Together with coordinates for the centre of the image plane and its
 * physical dimensions, this completely specifies the geometry of the scene.
 *
 * The actual implementation of these functions is in raycast.c.
 */

#include <math.h>
#include <stdlib.h>

#define Colour_size sizeof(float[3])
typedef float *Colour; // A Colour is an array of 3 floats: one per RGB channel.
                       // Because we need this to live on the heap, it has
                       // pointer type.

typedef double Vector[3]; // A point in space or a direction.

struct ImagePlane {
    struct {
        unsigned rows, cols;
    } resol;

    struct {
        double dims[2]; // Dimensions of the plane
        Vector centre;
        Vector tangent[2]; // Two tangent vectors. See note below about buffer.
    } geom;

    Colour **buff; // Buffer containing the contents of the image
    // The buffer should be accessed like:
    //          this.buff[row][col]
    // Note: that the image buffer is row-major: the "x" tangent vector points
    // down, in the direction that rows increase, and the "y" tangent vector
    // points right, in the direction that columns increase.
};

struct VoxelCube {
    struct {
        unsigned x, y, z;
    } resol;

    struct {
        // The geometry of the cube is trivial but we may as well store its
        // orientation for consistency and generalisability:
        double dims[3]; // dimensions of the cube
        Vector centre;
        Vector orient[3]; // Orientation of the cube
    } geom;

    Colour ***buff; // Buffer containing the contents of the voxel cloud
    // A voxel is just a colour.
    // The buffer should be accessed like:
    //          this.buff[row][col][layer]
};

// Construct and return a new unit cube:
struct VoxelCube new_unit_cube(unsigned x, unsigned y, unsigned z);

// Free the memory in the image buffer of the unit cube:
void free_unit_cube(struct VoxelCube unitcube);

// Construct and return a new image plane:
// We need to pass a reference cube, a distance from its centre, and two angles
// in polar coordinates. The function will return an image plane that is tangent
// to the bounding "sphere": i.e. oriented so that the vector connecting the
// cube's centre with the plane's centre is normal to the plane.
struct ImagePlane new_image_plane(struct VoxelCube ref_cube,
                                  double range, double theta, double phi,
                                  unsigned rows, unsigned cols);

// Free the memory in the image buffer of the image plane:
void free_image_plane(struct ImagePlane plane);

#endif //RAYCAST_H
