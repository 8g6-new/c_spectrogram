#ifndef SPECTRAL_FEATURES_H
#define SPECTRAL_FEATURES_H

    #include <stddef.h>
    #include <stdbool.h>
    #include <math.h>
    #include <fftw3.h>
    #include <string.h>
    #include "audio_io.h"

    typedef struct {
        float sample_rate;
        size_t output_size;
        size_t num_frequencies;
        float *frequencies;
        float *magnitudes;
        float *phases;
        float fft_times;
        float *phasers;
    } stft_d;



    stft_d stft(audio_data *audio,size_t window_size, size_t hop_size,float *window_values);
    void free_stft(stft_d *result);


    void calculate_frequencies(stft_d *result, size_t window_size, float sample_rate);
    static void init_fft_output(stft_d *result, unsigned int window_size, unsigned int hop_size, unsigned int num_samples);
    
    void window_function(float *window_values, size_t window_size, const char *window_type);
    void mel_filter(float min_f, float max_f, size_t n_filters, float sr, size_t fft_size, float *filter);

    static double hz_to_mel(double f, float mid_f);
    static double mel_to_hz(double m, float mid_f);
    void calculate_frequencies(stft_d *stft, size_t window_size, float sample_rate);
    size_t safe_diff(size_t a, size_t b);
#endif