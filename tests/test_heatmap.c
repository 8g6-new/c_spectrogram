#include "../src/png_tools/png_tools.h"


float random_float(float mean, float stddev) {
    float u1 = ((float)rand() / (RAND_MAX));
    float u2 = ((float)rand() / (RAND_MAX));
    float z = sqrtf(-2.0f * logf(u1)) * cosf(2.0f * M_PI * u2);
    return mean + stddev * z;
}




int main() {
    unsigned char bg_clr[4]  =  {0, 0, 50, 255};

    static const size_t w = 512, h = 512;

    heatmap_t *hm = heatmap_new(w, h);

    srand((unsigned int)time(NULL));

    for (unsigned i = 511; i>0; i--) {
        heatmap_add_weighted_point(hm, 0, i, i);
    }

    unsigned char *heatmap = (unsigned char *)malloc(w * h * 4);

    memset(heatmap, 0, w * h * 4);

    if (!heatmap) {
        perror("malloc");
        exit(1);
    }
    

    heatmap_render_default_to(hm, heatmap);
    add_bg(heatmap,w,h,bg_clr);
    heatmap_free(hm);

    

    unsigned char *resized_image = resize_image(heatmap, w, h, new_width, new_height);

    save_png("heatmap_resized.png", resized_image, new_width, new_height);

    // Free the memory
    free(heatmap);
    free(resized_image);

    return 0;
}