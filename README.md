raycaster
=========

This program is a simple test program that will eventually be used to make
animations using volumetric raycasting to visualise scientific data.

It compiles with

    cc pixel.c -c -o pixel.a
    cc raycast.c -c -o raycast.a
    cc pixel.a raycast.a render.c -o renderer

Running `./renderer` then produces a PPM image in binary format.
