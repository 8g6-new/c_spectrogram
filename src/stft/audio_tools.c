#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>


clock_t start, end;
float cpu_time_used;
#define EXP 6 // microseconds

float mag = 0.0, phase = 0.0;

#include "audio_tools.h"
#include "../libheatmap/heatmap.h"
#include "../png_tools/png_tools.h"

unsigned char* heatmap_get_vars(size_t w, size_t h, heatmap_t **hm_out) {
    unsigned char *heatmap = malloc(sizeof(unsigned char) * w * h * 4);
    if (!heatmap) {
        perror("malloc heatmap");
        return NULL;
    }
    
    memset(heatmap, 0, w * h * 4);

    heatmap_t *hm = heatmap_new(w,h);
    
    if (!hm) {
        fprintf(stderr, "Failed to create heatmap.\n");
        free(heatmap);
        return NULL;
    }

    *hm_out = hm;

    return heatmap;
}


inline void free_stft(stft_d *result) {
    free(result->phasers);
    result->phasers = NULL;

    free(result->magnitudes);
    result->magnitudes = NULL;

    free(result->phases);
    result->phases = NULL;
}

inline static double hz_to_mel(double f, float mid_f) {
    return mid_f * log10(1 + f / 700.0);
}

inline static double mel_to_hz(double m, float mid_f) {
    return 700.0 * (pow(10, m / mid_f) - 1);
}

inline static void init_fft_output(stft_d *result,unsigned int window_size,unsigned int hop_size,unsigned int num_samples) {
    result->phasers         = NULL;
    result->magnitudes      = NULL;
    result->phases          = NULL;
    
    result->fft_times       = 0.0;
    result->num_frequencies = (size_t) window_size / 2;
    result->output_size     = (size_t) (num_samples/ hop_size);

    result->magnitudes      = (float*) malloc(result->num_frequencies * result->output_size * sizeof(float));
    result->phases          = (float*) malloc(result->num_frequencies * result->output_size * sizeof(float));
    result->phasers         = (float*) malloc(result->num_frequencies * 2 * result->output_size * sizeof(float));

    if (!result->phasers  || !result->magnitudes || !result->phases) {
        fprintf(stderr, "Memory allocation failed\n");
        free(result->phasers);
        free(result->magnitudes);
        free(result->phases);
    }
}


inline audio_data read_file(const char *filename){

    audio_data audio = {0};
    
    SNDFILE *file;
    SF_INFO sf_info;

    file = sf_open(filename, SFM_READ, &sf_info);
    if (!file) {
        fprintf(stderr, "Error opening file\n");
        return audio ;
    }

    audio.num_samples = (size_t)sf_info.frames * sf_info.channels;
    audio.samples     = (float*)malloc(audio.num_samples * sizeof(float));

    if (!audio.samples) {
        fprintf(stderr, "Memory allocation failed\n");
        free(audio.samples);
        audio.samples = NULL;
        sf_close(file);
        return audio ;
    }

    if (sf_readf_float(file,audio.samples, sf_info.frames) < sf_info.frames) {
        fprintf(stderr, "Error reading audio data\n");
        free(audio.samples);
        audio.samples = NULL;
        sf_close(file);
        return audio ;
    }

    audio.sample_rate = sf_info.samplerate;
    audio.channels    = sf_info.channels;
    audio.frames      = sf_info.frames;

    sf_close(file);

    return audio;
}




inline size_t safe_diff(size_t a, size_t b) {
    return a==b ? 1 : a-b;
}

void mel_filter(float min_f, float max_f, size_t n_filters, float sr, size_t fft_size, float *filter) {
    int fft_d = (int)fft_size / 2;

    float mid = 2595.0;

    float mel_min = hz_to_mel(min_f, mid);
    float mel_max = hz_to_mel(max_f, mid);

    float mel_step = (mel_max - mel_min) / (n_filters + 1);

    size_t freq_bins[n_filters + 2];
    for (size_t i = 0; i < n_filters + 2; i++) {
        float hz = mel_to_hz(mel_min + (mel_step * i), mid);
        freq_bins[i] = (size_t)round(hz * fft_size / sr);
    }

    for (size_t m = 0; m < n_filters; m++) {
        size_t start_bin = freq_bins[m];
        size_t mid_bin = freq_bins[m + 1];
        size_t end_bin = freq_bins[m + 2];

        for (size_t f = start_bin; f <= mid_bin; f++) {
            filter[m * fft_d + f] = (float)(f - start_bin) / safe_diff(mid_bin,start_bin);
        }

        for (size_t f = mid_bin + 1; f <= end_bin; f++) {
            filter[m * fft_d + f] = (float)(end_bin - f) / safe_diff(end_bin,mid_bin);
        }
    }
}


void window_function(float *window_values, size_t window_size, const char *window_type) {
    if (strcmp(window_type, "hann") == 0) {
        for (size_t i = 0; i < window_size; i++) {
            window_values[i] = 0.5 * (1 - cos(2 * M_PI * i / (window_size - 1)));
        }
    } else if (strcmp(window_type, "hamming") == 0) {
        for (size_t i = 0; i < window_size; i++) {
            window_values[i] = 0.54 - 0.46 * cos(2 * M_PI * i / (window_size - 1));
        }
    } else if (strcmp(window_type, "blackman") == 0) {
        for (size_t i = 0; i < window_size; i++) {
            window_values[i] = 0.42 - 0.5 * cos(2 * M_PI * i / (window_size - 1)) + 0.08 * cos(4 * M_PI * i / (window_size - 1));
        }
    } else if (strcmp(window_type, "bartlett") == 0) {
        for (size_t i = 0; i < window_size; i++) {
            window_values[i] = 1.0 - fabs((i - (window_size - 1) / 2.0) / ((window_size - 1) / 2.0));
        }
    } else if (strcmp(window_type, "blackman-harris") == 0) {
        for (size_t i = 0; i < window_size; i++) {
            window_values[i] = 0.35875 - 0.48829 * cos(2 * M_PI * i / (window_size - 1)) + 
                               0.14128 * cos(4 * M_PI * i / (window_size - 1)) - 
                               0.01168 * cos(6 * M_PI * i / (window_size - 1));
        }
    } else if (strcmp(window_type, "flattop") == 0) {
        for (size_t i = 0; i < window_size; i++) {
            window_values[i] = 1 - 1.93 * cos(2 * M_PI * i / (window_size - 1)) + 
                                 1.29 * cos(4 * M_PI * i / (window_size - 1)) - 
                                 0.388 * cos(6 * M_PI * i / (window_size - 1)) + 
                                 0.032 * cos(8 * M_PI * i / (window_size - 1));
        }
    } else if (strcmp(window_type, "gaussian") == 0) {
        float sigma = 0.4;  // Standard deviation
        for (size_t i = 0; i < window_size; i++) {
            window_values[i] = exp(-0.5 * pow((i - (window_size - 1) / 2.0) / (sigma * (window_size - 1) / 2.0), 2));
        }
    } else if (strcmp(window_type, "kaiser") == 0) {
        float alpha = 3.0;  // Shape parameter
        float denominator = 1.0 / (tgamma(alpha + 1) * tgamma(1 - alpha));
        for (size_t i = 0; i < window_size; i++) {
            float ratio = (i - (window_size - 1) / 2.0) / ((window_size - 1) / 2.0);
            window_values[i] = tgamma(alpha + 1) * tgamma(1 - alpha) * denominator * exp(alpha * sqrt(1 - ratio * ratio));
        }
    } else {
        fprintf(stderr, "Unknown window type. Using rectangular window.\n");
        for (size_t i = 0; i < window_size; i++) {
            window_values[i] = 1.0;
        }
    }
}

void calculate_frequencies(float *frequencies, size_t window_size, float sample_rate) {
    size_t half_window_size = window_size / 2;
    float scale = sample_rate / window_size;
    for (size_t i = 0; i < half_window_size; i++) {
        frequencies[i] = i * scale;
    }
}

//fftwf : 32bit vs fftw 64bit vs fftwq:  128bit

inline stft_d stft(float *samples,size_t window_size, size_t hop_size,uint64_t num_samples,float *window_values) {
    
    stft_d result = {0};

    init_fft_output(&result,window_size,hop_size,num_samples);

    fftwf_complex *in  = (fftwf_complex*)fftwf_malloc(sizeof(fftwf_complex) * window_size);
    fftwf_complex *out = (fftwf_complex*)fftwf_malloc(sizeof(fftwf_complex) * window_size);
    fftwf_plan p       = fftwf_plan_dft_1d((int)window_size, in, out, FFTW_FORWARD, FFTW_ESTIMATE);


    for (size_t i = 0; i < result.output_size; i++) {
        clock_t start = clock();

        size_t start_idx = i * hop_size;
        if (start_idx + window_size > num_samples) break;

        for (size_t j = 0; j < window_size; j++) {
            in[j][0] = samples[start_idx + j] * window_values[j];
            in[j][1] = 0.0;
        }

        fftwf_execute(p);

        clock_t end = clock();

  
        for (size_t j = 0; j < result.num_frequencies; j++) {

            size_t index = i * (result.num_frequencies * 2) + (j * 2); // number of samples = 2 * f_max and j* 2 ( real , imag)

            result.phasers[index]     = out[j][0];
            result.phasers[index + 1] = out[j][1];

            float mag   = sqrt(out[j][0] * out[j][0] + out[j][1] * out[j][1]);
            float phase = atan2(out[j][1], out[j][0]);

            result.magnitudes[i * result.num_frequencies + j] = mag;
            result.phases[i * result.num_frequencies + j]     = phase;
        }

        result.fft_times            += ((float)(end - start)) * pow(10, EXP) / CLOCKS_PER_SEC;
    }

    fftwf_destroy_plan(p);
    p = NULL;
    fftwf_free(in);
    in = NULL;

    fftwf_free(out);
    out=NULL;

    free(window_values);
    window_values = NULL;

    return result;
}


inline void spectrogram(stft_d *result,const char *output_file,unsigned char bg_clr[4],bool db){
    size_t w = result->output_size;
    size_t h = result->num_frequencies;

    heatmap_t *hm = NULL;
    unsigned char *heatmap = heatmap_get_vars(w,h, &hm);

    for (size_t t = 0; t < w; t++) {
        for (size_t f = 0; f < h; f++) {
            size_t inverted_axis = (h-f-1); 
            float mag   = result->magnitudes[t*h*inverted_axis];
            heatmap_add_weighted_point(hm,t,f,(float)db*log10(mag*mag+1) + (float)!db * mag); //brachless
        }
    }
   
    heatmap_render_default_to(hm, heatmap);

    add_bg(heatmap, w,h, bg_clr);
    heatmap_free(hm);

    save_png(output_file,heatmap, w,h);
}

inline void mel_spectrogram(stft_d *result,const char *output_file,unsigned char bg_clr[4],bool db,size_t num_filters,float *mel_filter_bank,float *mel_vales){
    
    size_t w = result->output_size;
    size_t h = result->num_frequencies;

    heatmap_t *hm = NULL;
    unsigned char *heatmap = heatmap_get_vars(w,num_filters,&hm);

    for (size_t t = 0; t < w; t++) {
        for (size_t mel = 0; mel < num_filters; mel++){
            float sum = 0;  
            for (size_t f = 0; f < h; f++) {
                size_t inverted_axis = (h-f-1); 

                float mag   = result->magnitudes[t*h+inverted_axis];

                sum+=mel_filter_bank[(num_filters - mel - 1) * h+inverted_axis] * mag * mag;
            }
            sum = (float)db*sum + (float)!db * log10(sum+1);
            
            mel_vales[t*num_filters+(num_filters - mel - 1)] = sum;

            heatmap_add_weighted_point(hm, t, mel,sum); 
        }
    }

    heatmap_render_default_to(hm, heatmap);

    add_bg(heatmap, w,num_filters, bg_clr);
    heatmap_free(hm);

    save_png(output_file,heatmap, w,num_filters);
}


inline void MFCC(float *mel_values,const char *output_file,size_t w,unsigned char bg_clr[4],bool db,size_t num_filters,size_t num_coff){
    
    heatmap_t *hm = NULL;
    unsigned char *heatmap = heatmap_get_vars(w,num_coff,&hm);

    for (size_t t = 0; t < w; t++) {
        for(size_t n=0;n<num_coff;n++){
            float sum = 0;
            for (size_t mel = 0; mel < num_filters; mel++){   
                sum+=(cos((M_PI/num_filters)*(mel+1/2)*n))  * mel_values[t*num_filters + (num_filters - mel - 1)];
            }
            heatmap_add_weighted_point(hm, t,n,pow(sum*sum,2));
        }
    }

    heatmap_render_default_to(hm, heatmap);

    add_bg(heatmap, w,num_coff, bg_clr);
    heatmap_free(hm);

    save_png(output_file,heatmap, w,num_coff);
}