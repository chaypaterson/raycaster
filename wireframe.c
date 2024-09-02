#include "raycast.h"

// Wireframe geometry used for drawing e.g. axes:
// NOT CURRENTLY FUNCTIONAL

void project(Vector point, struct ImagePlane plane) {
    // move the vector "point" onto the image plane projectively
    double ndotc = dot(plane.geom.normal, plane.geom.centre);
    double ndotr = dot(plane.geom.normal, point);
    //double ndote = dot(plane.geom.eye, point); // unused

    for (char axis = 0; axis < 3; ++axis) {
        // parallel projection:
        point[axis] += (ndotc - ndotr) * plane.geom.normal[axis];
        // TODO with the eye at a finite distance this becomes:
        // += (ndotc - ndotr) * (eye[axis]-point[axis]) / (ndote - ndotr);
    }
}

char is_on_plane(Vector point, struct ImagePlane plane) {
    // check if a point will be visible on the plane
    char test = 1; // really a bool. TODO?
    // each component should be less than half the width of the axis:
    for (char axis = 0; axis < 2; ++axis) { // watch out, 2D!
        double compt = dot(plane.geom.tangent[axis], point);
        test = test && (compt > -plane.geom.dims[axis] * 0.5);
        test = test && (compt < +plane.geom.dims[axis] * 0.5);
    }

    return test;
}

// TODO a voxel-free "draw axes" method? just project axes onto the plane? and
// give them colours?
void draw_axes(struct ImagePlane plane) {
    // draw an x axis, y axis, and z axis in colours from the origin
    // call me after wiping the image buffer and before raycasting
    // each axis should be 1 unit long
    double tmax = 1.0;
    double dt = plane.geom.dims[0] * 1.0 / plane.resol.rows; // guess
    // Directions to draw in:
    Vector directions[3] = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
    // R G B values for each axis:
    float colours[3][3] = {{1.0, 0, 0}, {0, 1.0, 0}, {0, 0, 1.0}};

    for (char axis = 0; axis < 3; ++axis) {
        Vector point = {0, 0, 0};
        // The axes will lie along lines: so their projections on the image
        // plane should also be lines. We can project the point and the
        // directions now to ensure we will be working in the plane:

        project(point, plane);
        project(directions[axis], plane);
        for (double t = 0; t < tmax; t += dt) {
            for (char axis2 = 0; axis2 < 3; ++axis2) {
                point[axis2] = t * directions[axis][axis2];
            }
            // bounds checking:
            if (!is_on_plane(point, plane)) break;

            // get pixel coordinates:
            double alpha = dot(plane.geom.tangent[0], point);
            unsigned row = (alpha + 0.5) * plane.resol.rows;
            double beta = dot(plane.geom.tangent[1], point);
            unsigned col = (beta + 0.5) * plane.resol.cols;
            printf("%d %d\n", row, col);

            // draw colour:
            memcpy(plane.buff[row][col], colours[axis], Colour_size);
        }
    }
}


