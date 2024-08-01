#ifndef STFT_H
#define STFT_H

#include <sndfile.h>
#include <fftw3.h>

    typedef struct {
        double frequency;
        double magnitude;
    } FrequencyMagnitudePair;

    typedef struct {
        double *phasers;
        int output_size;
        double *frequencies;
        int num_frequencies;
    } STFTResult;

    typedef struct {
        float top_freq;
        float mean_freq;
        float weighted_mean;
        float std_freq_dev;
        double time;
        float sdfbpf;
        float sdfbmif;
        float sdfbmxf;
    } STFTStats;


    STFTStats stft_weighted_avg(double *output, int num_frequencies, double *frequencies,unsigned short int num_samples,unsigned short int peak,unsigned short int min,unsigned short int max);
    int compare_by_magnitude(const void *a, const void *b);
    void apply_window_function(double *window_values, int window_size, const char *window_type);
    void calculate_frequencies(double *frequencies, int window_size, double sample_rate);
    STFTResult stft(const char *filename, int window_size, int hop_size, const char *window_type,double fft_times[]);
    void mel_filterbank(int num_filters, int fft_size, float sample_rate, float **filterbank);

#endif