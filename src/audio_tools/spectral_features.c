
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


 /**
 * @brief Zeroth order modified Bessel function of the first kind (I₀)
 */
static double bessel_i0(double x) {
    double sum = 1.0;
    double y = x * x / 4.0;
    double term = y;
    int k = 1;
    while (term > 1e-10 * sum) {
        sum += term;
        k++;
        term *= y / (k * k);
    }
    return sum;
}
 


inline size_t hz_to_index(size_t num_freq, size_t sample_rate, float f) {
    return (size_t)((num_freq * f * 2) / sample_rate);
}

bool is_power_of_two(size_t x) {
    return x && ((x & (x - 1)) == 0);
}

void *aligned_alloc_batch(size_t batch, size_t count, bool zero_init) {
    if (!is_power_of_two(batch)) {
        ERROR("Alignment (%zu) must be a power of two.", batch);
        return NULL;
    }

    size_t total_size = ((count + batch - 1) / batch) * batch;
    void *ptr         = aligned_alloc(batch, total_size);

    if (!ptr) {
        ERROR("aligned_alloc failed for size %zu with alignment %zu.", total_size, batch);
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

    if (hop_size == 0)
        WARN("Hop size should be >= 0");

    if (window_size == 0)
        WARN("Window size should be >= 0");

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


double hz_to_mel(double hz)       { return 2595.0 * log10(1.0 + hz / 700.0); }
double mel_to_hz(double mel)      { return 700.0 * (pow(10.0, mel / 2595.0) - 1.0); }

double hz_to_bark(double hz)      { return 6.0 * asinh(hz / 600.0); }
double bark_to_hz(double bark)    { return 600.0 * sinh(bark / 6.0); }

double hz_to_erb(double hz)       { return 21.4 * log10(4.37e-3 * hz + 1.0); }
double erb_to_hz(double erb)      { return (pow(10.0, erb / 21.4) - 1.0) / 4.37e-3; }

double hz_to_chirp(double hz)     { return log2(hz + 1.0); }
double chirp_to_hz(double chirp)  { return pow(2.0, chirp) - 1.0; }

double hz_to_cam(double hz)       { return 45.5 * log10(1.0 + hz / 700.0); }
double cam_to_hz(double cam)      { return 700.0 * (pow(10.0, cam / 45.5) - 1.0); }

double hz_to_log10(double hz)     { return log10(hz + 1.0); }
double log10_to_hz(double val)    { return pow(10.0, val) - 1.0; }


void print_melbank(const filter_bank_t *v) {
    if (!v || !v->weights || v->size == 0) {
        printf("melbank is empty or NULL.\n");
        return;
    }

    printf("melbank size: %zu\n", v->size);

    for (size_t i = 0; i < v->size;i++) 
        printf("Index %zu -> x: %zu, Weight: %.6f\n",i, v->freq_indexs[i],v->weights[i]);
}


filter_bank_t gen_filterbank(filter_type_t type,
                             float min_f, float max_f,
                             size_t n_filters, float sr,
                             size_t fft_size, float *filter) {
    
    filter_bank_t non_zero = {
        .freq_indexs = NULL,
        .weights     = NULL,
        .size        = 0,
        .num_filters = n_filters
    };

    const size_t num_f     = fft_size / 2;
    const size_t avg_len   = n_filters * (fft_size * 3 / 512);
    non_zero.freq_indexs   = malloc(avg_len * sizeof(size_t));
    non_zero.weights       = malloc(avg_len * sizeof(float));
    
    double (*hz_to_scale)(double) = NULL;
    double (*scale_to_hz)(double) = NULL;

    switch (type) {
        case F_MEL:     hz_to_scale = hz_to_mel;     scale_to_hz = mel_to_hz; break;
        case F_BARK:    hz_to_scale = hz_to_bark;    scale_to_hz = bark_to_hz; break;
        case F_ERB:     hz_to_scale = hz_to_erb;     scale_to_hz = erb_to_hz; break;
        case F_CHIRP:   hz_to_scale = hz_to_chirp;   scale_to_hz = chirp_to_hz; break;
        case F_CAM:     hz_to_scale = hz_to_cam;     scale_to_hz = cam_to_hz; break;
        case F_LOG10:   hz_to_scale = hz_to_log10;   scale_to_hz = log10_to_hz; break;
        default:
            ERROR("Unknown filterbank type enum: %d", type);
            return non_zero;
    }

    double scale_min = hz_to_scale((double)min_f);
    double scale_max = hz_to_scale((double)max_f);
    double step      = (scale_max - scale_min) / (double)(n_filters + 1);

    double bin_edges[n_filters + 2];

    for (size_t i = 0; i < n_filters + 2; i++) {
        double scale = scale_min + step * (double)i;
        double hz    = scale_to_hz(scale);
        bin_edges[i] = hz * (double)(num_f) / (sr / 2.0);
    }

    for (size_t m = 1; m <= n_filters; m++) {
        double left   = bin_edges[m - 1];
        double center = bin_edges[m];
        double right  = bin_edges[m + 1];

        int k_start = (int)floor(left);
        int k_end   = (int)ceil(right);

        for (int k = k_start; k <= k_end; k++) {
            if (k < 0 || k >= (int)num_f) continue;

            double weight = 0.0;
            if (k < center) {
                weight = (k - left) / (center - left);
            } else if (k <= right) {
                weight = (right - k) / (right - center);
            } else {
                continue;
            }

            non_zero.weights[non_zero.size]     = (float)weight;
            filter[(m - 1) * num_f + k]         = (float)weight;
            non_zero.freq_indexs[non_zero.size] = k;
            non_zero.size++;
        }
    }

    non_zero.freq_indexs = realloc(non_zero.freq_indexs, non_zero.size * sizeof(size_t));
    non_zero.weights     = realloc(non_zero.weights,     non_zero.size * sizeof(float));

    return non_zero;
}


/**
 * @brief Generate window coefficients
 *
 * @param window_values Pointer to pre-allocated array of size `window_size`
 * @param window_size Length of the window
 * @param window_type Name of the window (e.g., "hann", "kaiser")
 */
void window_function(float *window_values, size_t window_size, const char *window_type) {
    const double N = (double)window_size;

    if (strcmp(window_type, "hann") == 0) {
        for (size_t i = 0; i < window_size; i++)
            window_values[i] = 0.5f * (1.0f - cosf(2.0f * M_PI * i / (N - 1.0f)));

    } else if (strcmp(window_type, "hamming") == 0) {
        for (size_t i = 0; i < window_size; i++)
            window_values[i] = 0.54f - 0.46f * cosf(2.0f * M_PI * i / (N - 1.0f));

    } else if (strcmp(window_type, "blackman") == 0) {
        for (size_t i = 0; i < window_size; i++) {
            double a0 = 0.42, a1 = 0.5, a2 = 0.08;
            double phase = 2.0 * M_PI * i / (N - 1.0);
            window_values[i] = a0 - a1 * cos(phase) + a2 * cos(2.0 * phase);
        }

    } else if (strcmp(window_type, "blackman-harris") == 0) {
        for (size_t i = 0; i < window_size; i++) {
            double phase = 2.0 * M_PI * i / (N - 1.0);
            window_values[i] = 0.35875
                             - 0.48829 * cos(phase)
                             + 0.14128 * cos(2.0 * phase)
                             - 0.01168 * cos(3.0 * phase);
        }

    } else if (strcmp(window_type, "bartlett") == 0) {
        for (size_t i = 0; i < window_size; i++)
            window_values[i] = 1.0f - fabsf((float)(i - (N - 1.0) / 2.0) / ((N - 1.0) / 2.0));

    } else if (strcmp(window_type, "flattop") == 0) {
        for (size_t i = 0; i < window_size; i++) {
            double phase = 2.0 * M_PI * i / (N - 1.0);
            window_values[i] = 1.0
                             - 1.93 * cos(phase)
                             + 1.29 * cos(2.0 * phase)
                             - 0.388 * cos(3.0 * phase)
                             + 0.028 * cos(4.0 * phase);
        }

    } else if (strcmp(window_type, "gaussian") == 0) {
        double sigma = 0.4;
        double denom = sigma * (N - 1.0) / 2.0;
        for (size_t i = 0; i < window_size; i++) {
            double x = (i - (N - 1.0) / 2.0) / denom;
            window_values[i] = exp(-0.5 * x * x);
        }

    } else if (strcmp(window_type, "kaiser") == 0) {
        double alpha = 3.0;  // Shape parameter, can be adjusted (typical: 3~8)
        double denom = bessel_i0(alpha);
        for (size_t i = 0; i < window_size; i++) {
            double ratio = 2.0 * i / (N - 1.0) - 1.0;
            window_values[i] = (float)(bessel_i0(alpha * sqrt(1.0 - ratio * ratio)) / denom);
        }

    } else {
    WARN("Unknown window type: %s. Using rectangular window.", window_type);
    for (size_t i = 0; i < window_size; i++)
        window_values[i] = 1.0f;
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

    long long elp = 0;
    long long temp = 0;

    for (size_t i = 0; i < output_size; i++) {
        const size_t start_idx = i * hop_size;
        for (size_t j = 0; j < window_size; j++) in[j][0] = mono[start_idx + j] * window_values[j];
        temp = get_time_us();
        fftwf_execute(plan);
        elp += (get_time_us() - temp);
        memcpy(&(result.phasers[i * window_size]), out,c_size);
    }

    double avg_elp  = (double)elp/output_size;

    printf("Avg time %f",avg_elp);


    FFT_bench(avg_elp,window_size);

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