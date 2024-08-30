raycaster
=========

This program is a simple test program that will eventually be used to make
animations using volumetric raycasting to visualise scientific data.

It compiles with

    cc pixel.c -c -o pixel.a
    cc raycast.c -c -o raycast.a
    cc pixel.a raycast.a render.c -o renderer

Running `./renderer` with no arguments then produces a series of PPM images in
binary format. This stack of PPM images can then be converted to a video with
e.g. FFmpeg.

The structure of render.c is deliberately simple and hackable, with most of the
work done by the pixel.c and raycast.c libraries. The API for these libraries is
detailed in pixel.h and raycast.h.

Note: compiling raycast.c with optimisation levels higher than -O1 under clang
results in unexpected behaviour and is not recommended. With GCC this doesn't
seem to happen, and -O3 is perfectly fine.

## Serialisation

The renderer can read voxels from files by invoking `./renderer --load
[yourfile.cube]`. A voxel cube file has the following format:

* File signature (first 6 bytes): "Voxel\n"
* Header (3 unsigned ints): resolution: rows, columns, and layers.
* The rest of the file (3 * rows * columns * layers * sizeof(float)): voxel colour data.

Currently voxel colour data is 3 channel floating point. This might change in
the future (alpha channels/opacity maybe) but is not a priority.

Once you have written your 3D colour data to a file, you can load that colour
data into the cube and render it by invoking `./renderer --load [yourfile.cube]`.
