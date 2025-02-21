#ifndef AUDIO_ANA_H
#define AUDIO_ANA_H

    #include <stdlib.h>
    #include <math.h>
    #include <sndfile.h>
    #include <fftw3.h>
    #include <stddef.h>
    #include <stdbool.h>
   
    #include "audio_io.h"


   
    typedef struct {
        float frequency;
        float magnitude;
    } freq_mag_pair;
    
     
    typedef struct {
        float *mean_mags;
        float *mean_freqs;
        float top_freq;
    } stft_mean;

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


    void calculate_frequencies(stft_d *result, size_t window_size, float sample_rate);

    static void init_fft_output(stft_d *result, unsigned int window_size, unsigned int hop_size, unsigned int num_samples);
    stft_d stft(audio_data *audio, size_t window_size, size_t hop_size, float *window_values);
    void free_stft(stft_d *result);


    stft_mean compute_stft_stats(stft_d *result,float f_min,float f_max,bool db);
    int compare_by_magnitude(const void *a, const void *b);



#endif