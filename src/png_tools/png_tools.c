#include "png_tools.h"


void save_png(const char *filename, const unsigned char *image, size_t width, size_t height) {
    FILE *fp = fopen(filename, "wb");
    if (!fp) {
        perror("fopen");
        exit(1);
    }

    png_structp png = png_create_write_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
    if (!png) {
        fprintf(stderr, "png_create_write_struct failed\n");
        fclose(fp);
        exit(1);
    }

    png_infop info = png_create_info_struct(png);
    if (!info) {
        fprintf(stderr, "png_create_info_struct failed\n");
        png_destroy_write_struct(&png, NULL);
        fclose(fp);
        exit(1);
    }

    if (setjmp(png_jmpbuf(png))) {
        fprintf(stderr, "Error during PNG creation\n");
        png_destroy_write_struct(&png, &info);
        fclose(fp);
        exit(1);
    }

    png_init_io(png, fp);
    png_set_IHDR(
        png,
        info,
        (png_uint_32)width, (png_uint_32)height,
        8,
        PNG_COLOR_TYPE_RGBA,
        PNG_INTERLACE_NONE,
        PNG_COMPRESSION_TYPE_DEFAULT,
        PNG_FILTER_TYPE_DEFAULT
    );

    png_write_info(png, info);

    png_bytep row = (png_bytep)malloc(4 * width);
    for (size_t y = 0; y < height; ++y) {
        memcpy(row, image + 4 * width * y, 4 * width);
        png_write_row(png, row);
    }

    printf("\n%s saved\n", filename);

    free(row);

    png_write_end(png, NULL);
    png_destroy_write_struct(&png, &info);
    fclose(fp);
}

unsigned char* resize_image(const unsigned char *original, size_t orig_width, size_t orig_height, size_t new_width, size_t new_height) {
    // Allocate memory for the resized image
    unsigned char *resized = (unsigned char *)malloc(new_width * new_height * 4);
    if (!resized) {
        perror("malloc");
        exit(1);
    }

    // Compute scaling factors
    float x_scale = (float)orig_width / new_width;
    float y_scale = (float)orig_height / new_height;

    // Pointers to optimize access
    const unsigned char *orig_row_ptr;
    unsigned char *resized_row_ptr = resized;

    // Size of a row in bytes for the original and resized images
    size_t orig_row_size = orig_width * 4;
    size_t resized_row_size = new_width * 4;

    // Resize the image
    for (size_t y = 0; y < new_height; ++y) {
        orig_row_ptr = original + (size_t)(y * y_scale) * orig_row_size;
        unsigned char *resized_pixel_ptr = resized_row_ptr;

        for (size_t x = 0; x < new_width; ++x) {
            size_t orig_x = (size_t)(x * x_scale) * 4;
            memcpy(resized_pixel_ptr, orig_row_ptr + orig_x, 4); // Copy pixel data
            resized_pixel_ptr += 4;
        }
        resized_row_ptr += resized_row_size;
    }
    
    printf("Image resized from (%ld x %ld) -> (%ld x %ld)", orig_width, orig_height,new_width,new_height);

    return resized;
}
