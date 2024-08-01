#include "../src/png_tools/png_tools.h"


float random_float(float mean, float stddev) {
    float u1 = ((float)rand() / (RAND_MAX));
    float u2 = ((float)rand() / (RAND_MAX));
    float z = sqrtf(-2.0f * logf(u1)) * cosf(2.0f * M_PI * u2);
    return mean + stddev * z;
}

int main() {
    static const size_t w = 512, h = 512, npoints = 1000;

    // Create the heatmap object with the given dimensions (in pixels).
    heatmap_t* hm = heatmap_new(w, h);

    // Seed the random number generator
    srand((unsigned int)time(NULL));

    // Add a bunch of random points to the heatmap.
    for (unsigned i = 0; i < npoints; ++i) {
        float x = random_float(0.5f * w, 0.5f / 3.0f * w);
        float y = random_float(0.5f * h, 0.25f * h);
        heatmap_add_point(hm, x, y);
    }

    // Allocate memory for the image
    unsigned char *image = (unsigned char *)malloc(w * h * 4);
    if (!image) {
        perror("malloc");
        exit(1);
    }

    // Initialize the image with black color
    memset(image, 255, w * h * 4);

    // Render the heatmap on top of the black background
    heatmap_render_default_to(hm, image);

    // Now that we've got a finished heatmap picture, we don't need the map anymore.
    heatmap_free(hm);

    // Resize the image to 128x128
    size_t new_width = 128;
    size_t new_height = 128;
    unsigned char *resized_image = resize_image(image, w, h, new_width, new_height);

    // Save the resized image as PNG
    save_png("heatmap_resized.png", resized_image, new_width, new_height);

    // Free the memory
    free(image);
    free(resized_image);

    return 0;
}
