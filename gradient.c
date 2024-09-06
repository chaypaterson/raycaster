#include "pixel.h"
#include <math.h>

// Print a palette to a 640x320 image

#define swatch1(theta) {\
COLOUR_MAX * (0.666 - 0.333 * cos((theta))),\
COLOUR_MAX * (0.666 + 0.167 * cos((theta)) - 0.288 * sin((theta))),\
COLOUR_MAX * (0.666 + 0.167 * cos((theta)) + 0.288 * sin((theta)))}

#define swatch2(theta) {\
COLOUR_MAX * (0.5 - 0.408 * cos((theta))),\
COLOUR_MAX * (0.5 + 0.204 * cos((theta)) - 0.353 * sin((theta))),\
COLOUR_MAX * (0.5 + 0.204 * cos((theta)) + 0.353 * sin((theta)))}

#define swatch3(theta) {\
COLOUR_MAX * (0.5 - 0.392 * cos((theta)) + 0.276 * sin((theta))),\
COLOUR_MAX * (0.5 + 0.087 * cos((theta)) - 0.122 * sin((theta))),\
COLOUR_MAX * (0.5 + 0.298 * cos((theta)) + 0.398 * sin((theta)))}

inline double round(double f) {
    // round to nearest 16th:
    return ((int)(f * 16)) * 1.0 / 16;
}

void gradient_test() {
    // Write out a test image:
    const char* imagefile = "testimg.ppm";
    Pixel BLACK = {0, 0, 0};
    // NB: COLOUR_MAX, 255, and 0xFF are synonymous on most systems.
    Pixel MAGENTA = {COLOUR_MAX, 0x00, COLOUR_MAX};
    Pixel CYAN = {0x00, COLOUR_MAX, COLOUR_MAX};
    Pixel RED = {COLOUR_MAX, 0x00, 0x00};

    int cols = 640;
    int rows = 320;

    // Open a file and write an image to it:

    FILE* img = fopen(imagefile, "w");

    putheader(img, cols, rows);

    // Draw the image:
    for (int row = 0; row < rows; ++row) {
        for (int col = 0; col < cols; ++col) {
            // Choose a colour:
            double f = col * 1.0/cols;
            double theta = 2 * M_PI * f;
            Pixel colour = swatch2(theta);

            putpixel(img, colour);
        }
    }

    fclose(img);
}

int main(void) {
    // gradient test:
    gradient_test();

    return 0;
}
