#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "src/stft/audio_tools.h"
#include "src/png_tools/png_tools.h"

#define EXP 6  // microseconds

float hm_freq = 0.0, std_mean = 0.0, std_dev = 0.0;

char find_exp(unsigned short int exp) {
    switch (exp) {
        case 3:
            return 'm';
        case 6:
            return 'u';
        default:
            return '\0';
    }
}

const char *th[] = {"st", "nd", "rd", "th"};

char check(int i, int lim) {
    return (i < lim - 1) ? '\n' : ' ';
}

// void free_memory(float **filterbank, size_t num_filters, unsigned char *heatmap, unsigned char *resized_image) {
//     if (filterbank) {
//         for (size_t i = 0; i < num_filters; ++i) {
//             free(filterbank[i]);
//         }
//         free(filterbank);
//     }
//     if (heatmap) free(heatmap);
//     if (resized_image) free(resized_image);
// }

int main(int argc, char *argv[]) {
    if (argc != 7) {
        fprintf(stderr, "Usage: %s <filename> <window_size> <hop_size> <window_type> <number_of_mel_banks>\n", argv[0]);
        return 1;
    }
    
    clock_t start = clock();

    const char *filename = argv[1];
    int window_size = atoi(argv[2]);
    int hop_size = atoi(argv[3]);
    const char *window_type = argv[4];
    unsigned short int num_filters = (unsigned short int)atoi(argv[5]);
    const char *output_file = argv[6];

    if (window_size <= 0 || hop_size <= 0 || num_filters <= 0) {
        fprintf(stderr, "Window size, hop size, and number of mel banks must be positive integers.\n");
        return 1;
    }

    // float **filterbank = (float **)malloc(num_filters * sizeof(float *));
    // if (!filterbank) {
    //     perror("malloc");
    //     return 1;
    // }

    // for (size_t i = 0; i < num_filters; ++i) {
    //     filterbank[i] = (float *)calloc(window_size, sizeof(float));
    //     if (!filterbank[i]) {
    //         perror("calloc");
    //         free_memory(filterbank, i, NULL, NULL);
    //         return 1;
    //     }
    // }

    // mel_filterbank(num_filters, window_size, 44000, filterbank);

    STFTResult result = stft(filename, window_size, hop_size, window_type);
    
    if (result.phasers == NULL || result.frequencies == NULL || result.magnitudes == NULL || result.phases == NULL || result.output_size == 0) {
        fprintf(stderr, "STFT computation failed. Check if the file exists and if memory allocation was successful.\n");
        // free_memory(filterbank, num_filters, NULL, NULL);
        return 1;
    }

    size_t w = result.output_size;
    size_t h = result.num_frequencies;

    heatmap_t *hm = heatmap_new(w, h);
    if (!hm) {
        fprintf(stderr, "Failed to create heatmap.\n");
        return 1;
    }

    unsigned char *heatmap = (unsigned char *)malloc(w * h * 4);
    if (!heatmap) {
        perror("malloc");
        heatmap_free(hm);
        return 1;
    }
    memset(heatmap, 0, w * h * 4);

    for (size_t i = 0; i < w; i++) {

        for (size_t j = 0; j < h; j++) {
            size_t index = i * h + (h - j - 1); 

            float phase = result.phases[index];
            float mag   = result.magnitudes[index];

            mag = (mag-result.mag[0])/result.mag[1];
            
            unsigned char n = (int)(((phase -  result.phase[0])) * 255 / (( result.phase[1]- result.phase[0])));

            // heatmap[i* j * 4 + 3] =  n; 
            heatmap_add_weighted_point(hm, i, j,mag);
        }
    }

    heatmap_render_default_to(hm, heatmap);

    unsigned char bg_clr[4] = {0, 0,0, 255};
    add_bg(heatmap, w, h, bg_clr);
    heatmap_free(hm);

    // size_t new_width = 256;
    // size_t new_height = 256;

    // // unsigned char *resized_image = resize_image(heatmap, w, h, new_width, new_height);

    // if (!resized_image) {
    //     perror("resize_image");
    //     return 1;
    // }
    
    save_png(output_file,heatmap, w,h);

    // free_memory(num_filters, heatmap,resized_image);

    clock_t end = clock();

    printf("%f ms\n",(float)(end - start) * pow(10,3) / CLOCKS_PER_SEC);

    free(heatmap);
    free(result.phasers);
    free(result.frequencies);
    free(result.magnitudes);
    free(result.phases);
    free(result.fft_times);

    return 0;
}
