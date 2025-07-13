
#include "../../headers/audio_tools/spectral_features.h"

/*
 * The MIT License (MIT)
 * 
 * Copyright © 2025 Devadut S Balan
 * 
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the “Software”), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 * 
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

 


inline size_t hz_to_index(size_t num_freq, size_t sample_rate, float f) {
    return (size_t)((num_freq * f * 2) / sample_rate);
}

bool is_power_of_two(size_t x) {
    return x && ((x & (x - 1)) == 0);
}

void *aligned_alloc_batch(size_t batch, size_t count, bool zero_init) {
    if (!is_power_of_two(batch)) {
        fprintf(stderr, "Error: alignment (%zu) must be a power of two.\n", batch);
        return NULL;
    }

    size_t total_size = ((count + batch - 1) / batch) * batch;
    void *ptr         = aligned_alloc(batch, total_size);

    if (!ptr) {
        fprintf(stderr, "Error: aligned_alloc failed for size %zu with alignment %zu.\n", total_size, batch);
        return NULL;
    }

    if (zero_init) {
        memset(ptr, 0, total_size);
    }

    return ptr;
}


void free_fft_plan(fft_d *fft) {
    if (fft) {
        if (fft->plan!=NULL) {
            fftwf_destroy_plan(fft->plan);
            fft->plan = NULL;
        }
        if (fft->in!=NULL) {
            fftwf_free(fft->in);
            fft->in = NULL;
        }
        if (fft->out!=NULL) {
            fftwf_free(fft->out);
            fft->out = NULL;
        }
    }
}



// void cleanup_fft_threads(fft_d *thread_ffts, const size_t num_threads) {
//     if (!thread_ffts) return;

//     for (size_t i = 0; i < num_threads; i++) {
//         free_fft_plan(&thread_ffts[i]);
//     }

//     free(thread_ffts);
//     thread_ffts = NULL;
// }

void free_stft(stft_d *result) {
    if (result) {
        free(result->phasers);
        free(result->magnitudes);
        free(result->frequencies);
        result->phasers = NULL;
        result->magnitudes = NULL;
        result->frequencies = NULL;
    }
}


bool init_fft_output(stft_d *result, unsigned int window_size, unsigned int hop_size, unsigned int num_samples) {

    if(hop_size==0)
    fprintf(stderr, "Hop size should be >= 0");

    if(window_size==0)
    fprintf(stderr, "window_size >= 0");

    if (!result) return false;

    result->phasers         = NULL;
    result->magnitudes      = NULL;
    result->frequencies     = NULL;


    result->num_frequencies = (size_t)(window_size / 2);
    result->output_size     = (size_t)(num_samples / hop_size);

    size_t magnitudes_size  = result->num_frequencies * result->output_size * sizeof(float);
    size_t phases_size      = result->num_frequencies * result->output_size * sizeof(float);
    size_t phasers_size     = result->num_frequencies * 2 * result->output_size * sizeof(float);

    result->magnitudes      = (float*)aligned_alloc_batch(32,magnitudes_size,false);
    result->phasers         = (float*)aligned_alloc_batch(32,phasers_size,false);

    return true;
}

inline  double hz_to_mel(double f, float mid_f) {
    return mid_f * log10(1 + f / 700.0);
}

inline  double mel_to_hz(double m, float mid_f) {
    return 700.0 * (pow(10, m / mid_f) - 1);
}

void print_melbank(const melbank_t *v) {
    if (!v || !v->weights || v->size == 0) {
        printf("melbank is empty or NULL.\n");
        return;
    }

    printf("melbank size: %zu\n", v->size);

    for (size_t i = 0; i < v->size;i++) 
        printf("Index %zu -> x: %zu, Weight: %.6f\n",i, v->freq_indexs[i],v->weights[i]);
}


melbank_t mel_filter(float min_f, float max_f, size_t n_filters, float sr, size_t fft_size, float *filter) {

    melbank_t non_zero      = { .freq_indexs = NULL, .weights = NULL, .size = 0 ,.num_filters=n_filters };
   
    const size_t avg_length = n_filters * (fft_size*3/512);
    non_zero.freq_indexs    = malloc( avg_length * sizeof(size_t));
    non_zero.weights        = malloc( avg_length * sizeof(float));


    
    size_t num_f            = (fft_size / 2);
    float mel_mid           = 2595.0f;
    float mel_min           = hz_to_mel(min_f, mel_mid);
    float mel_max           = hz_to_mel(max_f, mel_mid);
    float mel_step          = (mel_max - mel_min) / (n_filters + 1);
    
    float freq_bins[n_filters + 2];

    float mel,hz;
    
    for (size_t i = 0; i < n_filters + 2; i++) {
        mel          = mel_min + (mel_step * i);
        hz           = mel_to_hz(mel,mel_mid);
        freq_bins[i] = hz * (float)(num_f) / (sr / 2);
    }

    for (size_t m = 1; m <= n_filters; m++) {
        float left   = freq_bins[m - 1];
        float center = freq_bins[m];
        float right  = freq_bins[m + 1];
        
        int k_start = (int)floorf(left);
        int k_end   = (int)ceilf(right);

        for (int k = k_start; k <= k_end; k++) {  
            if (k < 0 || k >= num_f) continue;
            
            if (left <= k && k < center) {
                non_zero.weights[non_zero.size]     = (k - left) / (center - left);
                filter[(m - 1) * num_f + k]         =  non_zero.weights[non_zero.size];
                non_zero.freq_indexs[non_zero.size] = k;
                non_zero.size++;
            } else if (center <= k && k <= right) {
                non_zero.weights[non_zero.size]     = (right - k) / (right - center);
                filter[(m - 1) * num_f + k]         =  non_zero.weights[non_zero.size];
                non_zero.freq_indexs[non_zero.size] = k;
                non_zero.size++;
            }
        }
    }

    non_zero.freq_indexs = realloc(non_zero.freq_indexs, non_zero.size * sizeof(size_t));
    non_zero.weights     = realloc(non_zero.weights,     non_zero.size * sizeof(float));

    return non_zero;

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

void calculate_frequencies(stft_d *result, size_t window_size, float sample_rate) {
    size_t half_window_size = window_size / 2;
    float scale             = sample_rate / window_size;
    result->frequencies     = (float*) malloc(half_window_size * sizeof(float));
    for (size_t i = 0; i < half_window_size; i++) {
        result->frequencies[i] = i * scale;
    }
}

inline fft_d init_fftw_plan(const size_t window_size, const char *cache_dir) {
    fft_d fft = {0}; 

    char filename[256];
    snprintf(filename, sizeof(filename), "%s/%zu.wisdom", cache_dir, window_size);

    fft.in  = (fftwf_complex*)fftwf_malloc(sizeof(fftwf_complex) * window_size);
    fft.out = (fftwf_complex*)fftwf_malloc(sizeof(fftwf_complex) * window_size);


    if (!fft.in || !fft.out) {
        fprintf(stderr, "Memory allocation failed for FFTW buffers.\n");
        return fft;  
    }

    FILE *wisdom_file = fopen(filename, "r");
    if (wisdom_file) {
        if (fftwf_import_wisdom_from_file(wisdom_file)) {
            printf("Loaded optimized FFT plan: %s\n", filename);
            fft.plan = fftwf_plan_dft_1d((int)window_size, fft.in, fft.out, FFTW_FORWARD, FFTW_WISDOM_ONLY);
        } else {
            fprintf(stderr, "Error importing wisdom from file: %s\n", filename);
        }
        fclose(wisdom_file);
    }

    if (fft.plan == NULL) {
        printf("Cache not found or import failed. Creating FFT plan...\n");
        fft.plan = fftwf_plan_dft_1d((int)window_size, fft.in, fft.out, FFTW_FORWARD, FFTW_MEASURE);

        FILE *wisdom_out = fopen(filename, "w");
        if (wisdom_out) {
            fftwf_export_wisdom_to_file(wisdom_out);
            printf("Saved optimized FFT plan: %s\n", filename);
            fclose(wisdom_out);
        } else {
            fprintf(stderr, "Error saving wisdom to file: %s\n", filename);
        }
    }

    return fft;
}

stft_d stft(audio_data *audio, size_t window_size, size_t hop_size, float *window_values, fft_d *master_fft) {


    omp_set_num_threads(omp_get_max_threads());

    if (!audio || !window_values || !master_fft || !master_fft->plan) {
        stft_d error_result = {0};
        fprintf(stderr, "FFT Plan not found");
        return error_result;
    }

    stft_d result = {0};
    result.sample_rate = audio->sample_rate;


    init_fft_output(&result, window_size, hop_size, audio->num_samples);
    calculate_frequencies(&result, window_size, audio->sample_rate);

    audio->channels          = (audio->channels != 0) ? audio->channels : 1;
    const size_t channels    = audio->channels;
    const size_t length      = audio->num_samples/channels;
    const size_t output_size = (length - window_size) / hop_size + 1;
    result.output_size       = output_size;

    fftwf_complex *in       = master_fft->in;
    fftwf_complex *out      = master_fft->out;
    fftwf_plan plan         = master_fft->plan;

    float *mono              = malloc(length * sizeof(float)); 

    if (!mono) {
        stft_d error_result = {0};
        fprintf(stderr, "Memory allocation failed");
        return error_result;
    }

    const float   half   = 0.5f;
    const size_t  c_size = window_size * sizeof(fftwf_complex);

    if (channels == 1) memcpy(mono,audio->samples,length*sizeof(float));
    else if (channels == 2) {
        #pragma omp parallel for 
        for (size_t i = 0; i < length; i++) mono[i] = (audio->samples[i * 2] + audio->samples[i * 2 + 1]) * half;
    }

    memset(in, 0.0f, window_size * sizeof(fftwf_complex));

    for (size_t i = 0; i < output_size; i++) {
        const size_t start_idx = i * hop_size;
        for (size_t j = 0; j < window_size; j++) in[j][0] = mono[start_idx + j] * window_values[j];
        fftwf_execute(plan);
        memcpy(&(result.phasers[i * window_size]), out,c_size);
    }

    free(mono);

    const size_t total_size = output_size * result.num_frequencies;

    #pragma omp parallel for 
    for (size_t i = 0; i < total_size; i++) {
        const float real     = result.phasers[i * 2];
        const float imag     = result.phasers[i * 2 + 1];
        result.magnitudes[i] = sqrtf(real * real + imag * imag);
    }

    return result;
}