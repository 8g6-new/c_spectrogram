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
    int window_size = 4096;
    int hop_size = 512;
    const char *window_type = "hann";
    const char exp_symbol = find_exp(EXP);

    double fft_times[window_size];
    double real, imag, mag, phase;

    if (window_size <= 0 || hop_size <= 0) {
        fprintf(stderr, "Window size or hop size of samples must be positive integers.\n");
        return 1;
    }

    STFTResult result = stft(filename, window_size, hop_size, window_type, fft_times);

    if (result.phasers == NULL || result.frequencies == NULL) {
        fprintf(stderr, "STFT computation failed.\n");
        return 1;
    }

    for (int i = 0; i < result.output_size; i++) {
        printf("Frame %d\n", i + 1);

        for (int j = 0; j < result.num_frequencies; j++) {
            real = result.phasers[i * result.num_frequencies * 2 + j * 2];
            imag = result.phasers[i * result.num_frequencies * 2 + j * 2 + 1];
            mag = sqrt(real * real + imag * imag);
            phase = (mag > 0.0009) ? atan2(imag, real) : 0.0;

            printf("Phaser (%8.3f + %8.3f) Magnitude: %8.3f, Phase: %8.3f° (%8.3f Rads) @ %8.3f Hz\n",
                   real, imag, mag, phase * (180 / M_PI), phase, result.frequencies[j]);
        }

        printf("\nTime took to compute FFT of %d%s frame: %.3f %c%s\n", i + 1, th[i > 2 ? 3 : i], fft_times[i], exp_symbol, (exp_symbol != '\0') ? "s" : "");
        printf("\n");
    }

    free(result.phasers);
    free(result.frequencies);

    return 0;
}
