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

#define VChannels 3 // number of colour channels for voxels
#define Colour_size sizeof(float[VChannels])
typedef float *Colour; // A Colour is an array of 3 floats: one per RGB channel.
                       // Because we need this to live on the heap, it has
                       // pointer type instead of array type.

typedef double Vector[3]; // A point in space or a direction.

struct ImagePlane {
    // Resolution of image formed on plane:
    struct {
        unsigned rows, cols;
    } resol;

    // Geometry of the plane in the scene:
    struct {
        double dims[2]; // Dimensions of the plane
        Vector centre;
        Vector tangent[2]; // Two tangent vectors. See note below about buffer.
        Vector normal; // The vector normal to the plane
        Vector eye; // The eye from which rays are cast (optional)
    } geom;

    // "Film" is a buffer containing the contents of the image:
    Colour **film;
    // The film should be accessed like:
    //          this.film[row][col]
    // Note: that the image buffer is row-major: the "x" tangent vector points
    // down, in the direction that rows increase, and the "y" tangent vector
    // points right, in the direction that columns increase.
};

struct VoxelCube {
    // Resolution of the cube:
    struct {
        unsigned x, y, z; // TODO better names? how important?
    } resol;

    // The geometry of the cube: 
    struct {
        // This is currently trivial but we may as well store its
        // orientation for consistency and generalisability:
        double dims[3]; // dimensions of the cube
        Vector centre;
        Vector orient[3]; // Orientation of the cube
        double extinction; // Opacity of the cube. 1.0 = transparent,
        // 0.0=completely opaque. TODO really we should give individual voxels
        // an alpha channel.
        double circumradius; // Half the diagonal length of the cube
    } geom;

    Colour ***buff; // Buffer containing the contents of the voxel cloud
    // A voxel is a floating-point colour tuple.
    // The buffer should be accessed like:
    //          this.buff[row][col][layer]
};

// Construct and return a new unit cube with the specified resolution:
struct VoxelCube new_unit_cube(unsigned res_x, unsigned res_y, unsigned res_z);

// Free the memory in the image buffer of the unit cube:
void free_unit_cube(struct VoxelCube unitcube);

// Save and load cubes to and from files:
void save_cube(struct VoxelCube cube, char *filename);
struct VoxelCube load_cube(char *filename);

// Check the biggest channel value in the cube:
float maximum_colour_value(struct VoxelCube cube);

// Construct and return a new image plane with the specified resolution:
struct ImagePlane new_image_plane(unsigned rows, unsigned cols);

// Wipe the image plane, filling it with black pixels:
void wipe_plane(struct ImagePlane plane);

// Free the memory in the image buffer of the image plane:
void free_image_plane(struct ImagePlane plane);

// Quantise and write the image plane film to a file:
void save_image_plane(char* frameppm, struct ImagePlane plane, 
                      float exposure, float gamma);

// Reorient the image plane to view the reference cube from a given location in
// polar coordinates:
// We need to pass a reference cube, a distance from its centre, and two angles
// in polar coordinates. The function will return an image plane that is tangent
// to the bounding "sphere": i.e. oriented so that the vector connecting the
// cube's centre with the plane's centre is normal to the plane.
void orient_image_plane(struct ImagePlane* image_plane,
                        struct VoxelCube ref_cube,
                        double range, double theta, double phi);

// test if a point is inside a box:
_Bool is_inside_box(Vector point, struct VoxelCube box);

// shoot a ray along a vector and accumulate the colours of voxels inside a
// result:

void shoot_ray(Colour restrict result, Vector start, Vector dir,
               struct VoxelCube cube);

// cast rays from an image plane onto a voxel cube:

void raycast(struct ImagePlane image_plane, struct VoxelCube cube);

#endif //RAYCAST_H
