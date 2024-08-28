raycaster
=========

This program is a simple test program that will eventually be used to make
animations using volumetric raycasting to visualise scientific data.

It compiles with

    cc pixel.c -c -o pixel.a
    cc raycast.c -c -o raycast.a
    cc pixel.a raycast.a render.c -o renderer

Running `./renderer` then produces a series of PPM images in binary format.

Note: compiling raycast.c with optimisation levels higher than -O1 under clang
results in unexpected behaviour and is not recommended. With GCC this doesn't
seem to happen.

The renderer can read voxels from files. A voxel cube file has the following
format:

* File signature (6 bytes): "Voxel\n"
* Header (3 unsigned ints): resolution: rows, columns, and layers.
* The rest (3 * rows * columns * layers * sizeof(float)): voxel colour data.

Currently voxel colour data is 3 channel floating point. This might change in
the future (alpha channels/opacity maybe) but is not a priority.
