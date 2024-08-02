#ifndef STFT_H
#define STFT_H

#include <sndfile.h>
#include <fftw3.h>
#include <stddef.h>

    typedef struct {
        float frequency;
        float magnitude;
    } FrequencyMagnitudePair;

    typedef struct {
        float *phasers;
        size_t output_size;
        float *frequencies;
        unsigned int num_frequencies;
        float *magnitudes;
        float *phases;
        float *infos;
        unsigned short int *info_indexes;
        float *fft_times;
    } STFTResult;

    typedef struct {
        float top_freq;
        float mean_freq;
        float weighted_mean;
        float std_freq_dev;
        float time;
        float sdfbpf;
        float sdfbmif;
        float sdfbmxf;
    } STFTStats;


    STFTStats stft_weighted_avg(float*output, size_t num_frequencies, float*frequencies, unsigned short int num_samples, unsigned short int peak, unsigned short int min, unsigned short int max);
    int compare_by_magnitude(const void *a, const void *b);
    void apply_window_function(float *window_values, size_t window_size, const char *window_type);
    void calculate_frequencies(float *frequencies, size_t window_size, float sample_rate);
    STFTResult stft(const char *filename, size_t window_size, size_t hop_size, const char *window_type);
    void mel_filterbank(size_t num_filters, size_t fft_size, float sample_rate, float **filterbank);
#endif