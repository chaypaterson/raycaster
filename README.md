raycaster
=========

This program is a simple test program that will eventually be used to make
animations using volumetric raycasting to visualise scientific data.

## Building

It compiles (`-O3` optional but recommended) with

    cc [-O3] pixel.c -c -o pixel.a
    cc [-O3] raycast.c -c -o raycast.a
    cc [-O3] pixel.a raycast.a render.c -o renderer -lm

Running `./renderer` with no arguments then produces a series of PPM images in
binary format. This stack of PPM images can then be converted to a video with
e.g. FFmpeg.

The structure of render.c is deliberately simple and hackable, with most of the
work done by the pixel.c and raycast.c libraries. The API for these libraries is
detailed in pixel.h and raycast.h.

## Serialisation

The renderer can read voxels from files by invoking `./renderer --load
[yourfile.cube]`. A voxel cube file has the following format:

* File signature (first 6 bytes): "Voxel\n"
* Header (3 unsigned ints): resolution: rows, columns, and layers.
* The rest of the file (3 * rows * columns * layers * sizeof(float)): voxel colour data.

Currently voxel colour data is 3 channel floating point. This might change in
the future (alpha channels/opacity maybe) but is not a priority. Conventionally,
I put the file extension as `.cube` but this is arbitrary as long as the format
is as described above.

Once you have written your 3D colour data to a file, you can load that colour
data into the cube and render it by invoking `./renderer --load [yourfile.cube]`.

## Running

Run either `./renderer` for a demo mode, which will save an example voxel file
to disk, or `./renderer --load [yourfile.cube]` to load and make a video of that
file's contents.

NaN values will be rendered as black voxels or voxels with unexpected colours.
To keep this library simple, input is assumed to be nice and not contain NaNs.

## TODO

Currently, the image buffer and voxel buffers defined in raycast.h have type
`Colour **` and `Colour ***buff` -- it would be better to have types that
were contiguous in memory but could still be accessed as variable length arrays.

- The voxel buffer being contiguous in memory with a getter does not effect
  performance but the syntax `get_voxel(cube, row, col, lyr)` is not as clear as
  the syntax `cube.block[row][col][lyr]`.
