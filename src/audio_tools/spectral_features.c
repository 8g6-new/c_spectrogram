#include "../../headers/audio_tools/spectral_features.h"

#define CALCULATE_MAGNITUDE 1   
#define CALCULATE_PHASE     1  

void compute_magnitudes_and_phases_scalar(fftwf_complex *out, float *magnitudes, float *phases, const size_t offset, const size_t count) {
    for (size_t j = 0; j < count; j++) {
        const float real = out[j][0];
        const float imag = out[j][1];

        #ifdef CALCULATE_MAGNITUDE
            magnitudes[offset + j] = sqrtf(real * real + imag * imag);
        #endif

        #ifdef CALCULATE_PHASE
            phases[offset + j] = atan2f(imag, real);
        #endif
    }
}


inline  double hz_to_mel(double f, float mid_f) {
    return mid_f * log10(1 + f / 700.0);
}

inline  double mel_to_hz(double m, float mid_f) {
    return 700.0 * (pow(10, m / mid_f) - 1);
}

inline void free_fft(fft_d *fft) {
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

void cleanup_fft_threads(fft_d *thread_ffts, const size_t num_threads) {
    if (!thread_ffts) return;

    for (size_t i = 0; i < num_threads; i++) {
        free_fft(&thread_ffts[i]);
    }

    free(thread_ffts);
    thread_ffts = NULL;
}

void free_stft(stft_d *result) {
    if (result) {
        free(result->phasers);
        free(result->magnitudes);
        free(result->phases);
        free(result->frequencies);
        result->phasers = NULL;
        result->magnitudes = NULL;
        result->phases = NULL;
        result->frequencies = NULL;
    }
}


inline bool init_fft_output(stft_d *result, unsigned int window_size, unsigned int hop_size, unsigned int num_samples) {
    if (!result) return false;

    result->phasers         = NULL;
    result->magnitudes      = NULL;
    result->phases          = NULL;
    result->frequencies     = NULL;

    result->fft_times       = 0.0;
    result->num_frequencies = (size_t)(window_size / 2);
    result->output_size     = (size_t)(num_samples / hop_size);

    size_t magnitudes_size  = result->num_frequencies * result->output_size * sizeof(float);
    size_t phases_size      = result->num_frequencies * result->output_size * sizeof(float);
    size_t phasers_size     = result->num_frequencies * 2 * result->output_size * sizeof(float);

    result->magnitudes      = (float*)aligned_alloc(32, magnitudes_size);
    result->phases          = (float*)aligned_alloc(32, phases_size);
    result->phasers         = (float*)aligned_alloc(32, phasers_size);


    return true;
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


//fftwf : 32bit vs fftw 64bit vs fftwq:  128bit
inline stft_d stft(audio_data *audio, size_t window_size, size_t hop_size, float *window_values, fft_d *master_fft) {
    if (!audio || !window_values || !master_fft || !master_fft->plan) {
        stft_d error_result = {0};
        return error_result;
    }

    stft_d result                = {0};
    result.sample_rate           = audio->sample_rate;

    init_fft_output(&result, window_size, hop_size, audio->num_samples);
    calculate_frequencies(&result, window_size, audio->sample_rate);
    
    audio->channels              = (audio->channels != 0) ? audio->channels : 1;
    result.output_size          /= audio->channels;
    
    const size_t channels        = audio->channels;
    const size_t output_size     = result.output_size;
    const size_t num_frequencies = result.num_frequencies;
    const size_t num_samples     = audio->num_samples;

    size_t max_i = (num_samples >= window_size * channels) ? 
                    (num_samples - window_size * channels) / (hop_size ) : 0;
    max_i = (max_i < output_size) ? max_i : output_size;

    result.output_size           = max_i;

    const size_t num_threads     = omp_get_max_threads();
    fft_d *thread_ffts           = (fft_d*)malloc(num_threads * sizeof(fft_d));
    
    if (!thread_ffts) {
        fprintf(stderr, "Failed to allocate memory for thread FFT plans\n");
        return result;
    }
    
    #pragma omp parallel
    {
        int thread_id = omp_get_thread_num();
        
        thread_ffts[thread_id].in  = (fftwf_complex*)fftwf_malloc(window_size * sizeof(fftwf_complex));
        thread_ffts[thread_id].out = (fftwf_complex*)fftwf_malloc(window_size * sizeof(fftwf_complex));
        
        #pragma omp critical
        {
            thread_ffts[thread_id].plan = fftwf_plan_dft_1d(
                (int)window_size, 
                thread_ffts[thread_id].in, 
                thread_ffts[thread_id].out, 
                FFTW_FORWARD, 
                FFTW_WISDOM_ONLY 
            );
            
            if (!thread_ffts[thread_id].plan) {
                thread_ffts[thread_id].plan = fftwf_plan_dft_1d(
                    (int)window_size, 
                    thread_ffts[thread_id].in, 
                    thread_ffts[thread_id].out, 
                    FFTW_FORWARD, 
                    FFTW_ESTIMATE  
                );
            }
        }
    }

    
    #pragma omp parallel
    {
        int thread_id     = omp_get_thread_num();
        fft_d *thread_fft = &thread_ffts[thread_id];

        if (thread_fft->plan) {
            #pragma omp for 
            for (size_t i = 0; i < max_i; i++) {
                const size_t start_idx = i * hop_size ;

                for (size_t j = 0; j < window_size; j++) {
                    float sum = 0.0f;
                    for (size_t ch = 0; ch < channels; ch++) {
                        sum += audio->samples[start_idx + (j * channels) + ch];
                    }
                    thread_fft->in[j][0] = (sum * window_values[j]) / channels;
                    thread_fft->in[j][1] = 0.0f;
                }

                fftwf_execute(thread_fft->plan);

                const size_t offset = i * num_frequencies;
                compute_magnitudes_and_phases_scalar(thread_fft->out,result.magnitudes,result.phases,offset,num_frequencies);
            }
        }
    }

    cleanup_fft_threads(thread_ffts,num_threads);
    
    return result;
}