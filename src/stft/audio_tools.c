#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>

clock_t start, end;
double cpu_time_used;
#define EXP 6 // microseconds

double max = 0.0;
double real=0.0,imag=0.0,mag=0.0;

#include "audio_tools.h"

int compare_by_magnitude(const void *a, const void *b) {
    const FrequencyMagnitudePair *pair1 = (const FrequencyMagnitudePair *)a;
    const FrequencyMagnitudePair *pair2 = (const FrequencyMagnitudePair *)b;
    return (pair2->magnitude > pair1->magnitude) - (pair2->magnitude < pair1->magnitude);
}

STFTStats stft_weighted_avg(double *output, int num_frequencies, double *frequencies, unsigned short int num_samples,unsigned short int peak,unsigned short int min,unsigned short int max) {
    
    STFTStats stats;

    start = clock();

    if (num_frequencies <= 0) {
        fprintf(stderr, "Invalid number of frequencies %d\n",num_frequencies);
    }

    if (num_samples <= 0) {
        fprintf(stderr, "Invalid number of samples %d\n",num_samples);
    }

    FrequencyMagnitudePair *frequency_magnitude_pairs = (FrequencyMagnitudePair *)malloc(num_frequencies * sizeof(FrequencyMagnitudePair));

    if (!frequency_magnitude_pairs) {
        fprintf(stderr, "Memory allocation failed.\n");
    }

    
    for (int j = 0; j < num_frequencies; j++) {
        real = output[j * 2];
        imag = output[j * 2 + 1];
        mag  = sqrt(real * real + imag * imag);

        if (mag > max) max = mag;

        frequency_magnitude_pairs[j].frequency = frequencies[j];
        frequency_magnitude_pairs[j].magnitude = mag;
    }


    qsort(frequency_magnitude_pairs, num_frequencies, sizeof(FrequencyMagnitudePair), compare_by_magnitude);

    
    stats.weighted_mean  = 0.0;
    stats.std_freq_dev   = 0.0;
    stats.sdfbpf         = 0.0;
    stats.sdfbmif        = 0.0;
    stats.sdfbmxf        = 0.0;

    stats.top_freq = frequency_magnitude_pairs[0].frequency;

    if (num_samples > 0) {
         for (unsigned short int j = 0; j < num_samples; j++) {
            stats.mean_freq     += frequency_magnitude_pairs[j].frequency;
            stats.weighted_mean += frequency_magnitude_pairs[j].frequency * (frequency_magnitude_pairs[j].magnitude/max);
        }

         stats.mean_freq        = stats.mean_freq/num_samples;
        //  stats.weighted_mean    = stats.weighted_mean/num_samples;

        for (unsigned short int j = 0; j < num_samples; j++) {
            stats.std_freq_dev  += pow(frequency_magnitude_pairs[j].frequency - stats.mean_freq,2);
            stats.sdfbpf        += pow(frequency_magnitude_pairs[j].frequency - peak,2);
            stats.sdfbmif       += pow(frequency_magnitude_pairs[j].frequency - max,2);
            stats.sdfbmxf       += pow(frequency_magnitude_pairs[j].frequency - min,2);
        }

        stats.std_freq_dev     = sqrt(stats.std_freq_dev/num_samples);
        stats.sdfbpf           = sqrt(stats.sdfbpf/num_samples);
        stats.sdfbmif          = sqrt(stats.sdfbmif/num_samples);
        stats.sdfbmxf          = sqrt(stats.sdfbmxf/num_samples);

    }
    else{
       fprintf(stderr, "Number of samples is zero.\n %d",num_samples);
    }

    end        = clock();
    stats.time = ((double)(end - start))*pow(10,EXP) / CLOCKS_PER_SEC; 

    free(frequency_magnitude_pairs);

    return stats;
}

void mel_filterbank(int num_filters, int fft_size, float sample_rate, float **filterbank) {
    float f_min = 0;
    float f_max = sample_rate / 2;
    float center = 2595.0;
    float log_center = log10(1 + f_min / 700.0);
    float mel_min = center * log_center;
    float mel_max = center * log10(1 + f_max / 700.0);
    
    float *mel_points = (float *)malloc((num_filters + 2) * sizeof(float));
    float mel_step = (mel_max - mel_min) / (num_filters + 1);
    for (int i = 0; i < num_filters + 2; ++i) {
        mel_points[i] = mel_min + i * mel_step;
    }
    
    float *freq_points = (float *)malloc((num_filters + 2) * sizeof(float));
    int *bin_points = (int *)malloc((num_filters + 2) * sizeof(int));
    for (int i = 0; i < num_filters + 2; ++i) {
        freq_points[i] = 700 * (pow(10, mel_points[i] / center) - 1);
        bin_points[i] = round(freq_points[i] / (sample_rate / fft_size));
    }
    
    for (int m = 0; m < num_filters; ++m) {
        int start_bin = bin_points[m];
        int mid_bin = bin_points[m + 1];
        int end_bin = bin_points[m + 2];
        
        for (int k = start_bin; k < mid_bin; ++k) {
            filterbank[m][k] = (float)(k - start_bin) / (mid_bin - start_bin);
        }
        
        for (int k = mid_bin; k < end_bin; ++k) {
            filterbank[m][k] = (float)(end_bin - k) / (end_bin - mid_bin);
        }
    }
    
    // Free allocated memory
    free(mel_points);
    free(freq_points);
    free(bin_points);
}

void apply_window_function(double *window_values, int window_size, const char *window_type) {
    if (strcmp(window_type, "hann") == 0) {
        for (int i = 0; i < window_size; i++) {
            window_values[i] = 0.5 * (1 - cos(2 * M_PI * i / (window_size - 1)));
        }
    } else if (strcmp(window_type, "hamming") == 0) {
        for (int i = 0; i < window_size; i++) {
            window_values[i] = 0.54 - 0.46 * cos(2 * M_PI * i / (window_size - 1));
        }
    } else if (strcmp(window_type, "blackman") == 0) {
        for (int i = 0; i < window_size; i++) {
            window_values[i] = 0.42 - 0.5 * cos(2 * M_PI * i / (window_size - 1)) + 0.08 * cos(4 * M_PI * i / (window_size - 1));
        }
    } else if (strcmp(window_type, "bartlett") == 0) {
        for (int i = 0; i < window_size; i++) {
            window_values[i] = 1.0 - fabs((i - (window_size - 1) / 2.0) / ((window_size - 1) / 2.0));
        }
    } else if (strcmp(window_type, "blackman-harris") == 0) {
        for (int i = 0; i < window_size; i++) {
            window_values[i] = 0.35875 - 0.48829 * cos(2 * M_PI * i / (window_size - 1)) + 
                               0.14128 * cos(4 * M_PI * i / (window_size - 1)) - 
                               0.01168 * cos(6 * M_PI * i / (window_size - 1));
        }
    } else if (strcmp(window_type, "flattop") == 0) {
        for (int i = 0; i < window_size; i++) {
            window_values[i] = 1 - 1.93 * cos(2 * M_PI * i / (window_size - 1)) + 
                                 1.29 * cos(4 * M_PI * i / (window_size - 1)) - 
                                 0.388 * cos(6 * M_PI * i / (window_size - 1)) + 
                                 0.032 * cos(8 * M_PI * i / (window_size - 1));
        }
    } else if (strcmp(window_type, "gaussian") == 0) {
        double sigma = 0.4;  // Standard deviation
        for (int i = 0; i < window_size; i++) {
            window_values[i] = exp(-0.5 * pow((i - (window_size - 1) / 2.0) / (sigma * (window_size - 1) / 2.0), 2));
        }
    } else if (strcmp(window_type, "kaiser") == 0) {
        double alpha = 3.0;  // Shape parameter
        double denominator = 1.0 / (tgamma(alpha + 1) * tgamma(1 - alpha));
        for (int i = 0; i < window_size; i++) {
            double ratio = (i - (window_size - 1) / 2.0) / ((window_size - 1) / 2.0);
            window_values[i] = tgamma(alpha + 1) * tgamma(1 - alpha) * denominator * exp(alpha * sqrt(1 - ratio * ratio));
        }
    } else {
        fprintf(stderr, "Unknown window type. Using rectangular window.\n");
        for (int i = 0; i < window_size; i++) {
            window_values[i] = 1.0;
        }
    }
}

void calculate_frequencies(double *frequencies, int window_size, double sample_rate) {
    int half_window_size = window_size / 2;
    double scale = sample_rate / window_size;
    for (int i = 0; i < half_window_size; i++) {
        frequencies[i] = i * scale;
    }
}


STFTResult stft(const char *filename, int window_size, int hop_size, const char *window_type,double fft_times[]) {

    STFTResult result;

    result.phasers         = NULL;
    result.frequencies     = NULL;
    result.output_size     = 0;
    result.num_frequencies = window_size / 2;

    SNDFILE *file;
    SF_INFO sf_info;
    float *audio_data;
    int num_samples;

    file = sf_open(filename, SFM_READ, &sf_info);
    if (!file) {
        fprintf(stderr, "Error opening file\n");
        return result;
    }

    num_samples = sf_info.frames * sf_info.channels;
    audio_data = (float*)malloc(num_samples * sizeof(float));

    if (!audio_data) {
        fprintf(stderr, "Memory allocation failed\n");
        sf_close(file);
        return result;
    }

    if (sf_readf_float(file, audio_data, sf_info.frames) < sf_info.frames) {
        fprintf(stderr, "Error reading audio data\n");
        free(audio_data);
        sf_close(file);
        return result;
    }

    sf_close(file);

    fftw_complex *in  = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * window_size);
    fftw_complex *out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * window_size);
    fftw_plan p       = fftw_plan_dft_1d(window_size, in, out, FFTW_FORWARD, FFTW_ESTIMATE);

    double *window_values = (double*)malloc(window_size * sizeof(double));
    apply_window_function(window_values, window_size, window_type);

    result.frequencies = (double*)malloc(result.num_frequencies * sizeof(double));
    calculate_frequencies(result.frequencies, window_size, sf_info.samplerate);

    result.output_size = (num_samples - window_size) / hop_size + 1;
    result.phasers = (double*)malloc(result.num_frequencies * 2 * result.output_size * sizeof(double));

    for (int i = 0; i < result.output_size; i++) {
        
        start = clock();

        int start_idx = i * hop_size;
        if (start_idx + window_size > num_samples) break;
        

        for (int j = 0; j < window_size; j++) {
            in[j][0] = audio_data[start_idx + j] * window_values[j]; 
            in[j][1] = 0.0; 
        }

        fftw_execute(p);

        end = clock();


        for (int j = 0; j < result.num_frequencies; j++) {
            result.phasers[i * result.num_frequencies * 2 + j * 2    ] = out[j][0]; 
            result.phasers[i * result.num_frequencies * 2 + j * 2 + 1] = out[j][1]; 
        }

        fft_times[i] = ((double)(end - start))*pow(10,EXP) / CLOCKS_PER_SEC; 
    }

    fftw_destroy_plan(p);
    fftw_free(in);
    fftw_free(out);
    free(audio_data);
    free(window_values);

    return result;
}
