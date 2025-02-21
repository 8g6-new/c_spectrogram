#ifndef PNG_TOOLS_H
#define PNG_TOOLS_H

    #include <stddef.h>
    #include <stdlib.h>
    #include <stdio.h>
    #include <png.h>
    #include "../libheatmap/heatmap.h"
    #include <math.h>
    #include <string.h>

    void add_bg(unsigned char *image,  size_t width, size_t height, unsigned char color[4]);
    void save_png(const char *filename, const unsigned char *image, size_t width, size_t height);
    unsigned char* resize_image(const unsigned char *original, size_t orig_width, size_t orig_height, size_t new_width, size_t new_height);

#endif
