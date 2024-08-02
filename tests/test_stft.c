#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../src/stft/audio_tools.h"

#define EXP 6  // microseconds

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

int main() {
    const char *filename = "./tests/files/wind.wav";
    int window_size = 128;
    int hop_size = 1024;
    const char *window_type = "hann";
    const char exp_symbol = find_exp(EXP);

    if (window_size <= 0 || hop_size <= 0) {
        fprintf(stderr, "Window size or hop size of samples must be positive integers.\n");
        return 1;
    }

    STFTResult result = stft(filename, window_size, hop_size, window_type);

    if (result.phasers == NULL || result.frequencies == NULL || result.output_size == 0) {
        fprintf(stderr, "STFT computation failed. Check if the file exists and if memory allocation was successful.\n");
        return 1;
    }

    printf("Number of frames: %d\n", (int)result.output_size);
    printf("Number of frequencies: %d\n", (int)result.num_frequencies);

    for (size_t i = 0; i < result.output_size; i++) {
        printf("Frame %d\n", (int)i + 1);
        
        unsigned short int *indx = &result.info_indexes[i * 4];
        printf("\nMin mag @ %f Hz , Max mag @ %f\nMin phase @ %f,Max Phase @ %f\n\n",
            result.frequencies[indx[0]],  // Min magnitude index frequency
            result.frequencies[indx[1]],  // Max magnitude index frequency
            result.frequencies[indx[2]],  // Min phase index frequency
            result.frequencies[indx[3]]); // Max phase index frequency

        for(size_t j = 0; j < result.num_frequencies; j++) {
            size_t index = i * result.num_frequencies + j;
            if (index >= result.num_frequencies * result.output_size) {
                fprintf(stderr, "Index out of bounds. Check array sizes.\n");
                return 1;
            }

            float mag = result.magnitudes[index];
            float phase = result.phases[index];

            printf("Magnitude: %8.3f, Phase: %8.3f° (%8.3f Rads) @ %8.3f Hz\n",
                   mag, phase * (180 / M_PI), phase, result.frequencies[j]);
        }

        printf("\nTime took to compute FFT of %d%s frame: %.3f %c%s\n", (int)i + 1, th[i > 2 ? 3 : i], result.fft_times[i], exp_symbol, (exp_symbol != '\0') ? "s" : "");
        printf("\n");
    }

    free(result.phasers);
    free(result.frequencies);
    free(result.magnitudes);
    free(result.phases);
    free(result.fft_times);

    return 0;
}
