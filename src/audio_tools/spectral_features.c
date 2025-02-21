#include "../../headers/audio_tools/spectral_features.h"

inline static double hz_to_mel(double f, float mid_f) {
    return mid_f * log10(1 + f / 700.0);
}

inline static double mel_to_hz(double m, float mid_f) {
    return 700.0 * (pow(10, m / mid_f) - 1);
}



inline void free_stft(stft_d *result) {
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

inline void init_fft_output(stft_d *result, unsigned int window_size, unsigned int hop_size, unsigned int num_samples) {
    if (!result) return;

    result->phasers         = NULL;
    result->magnitudes      = NULL;
    result->phases          = NULL;
    result->frequencies     = NULL;
    
    result->fft_times       = 0.0;
    result->num_frequencies = (size_t)window_size / 2;
    result->output_size     = (size_t)(num_samples / hop_size);

    result->magnitudes      = (float*)malloc(result->num_frequencies * result->output_size * sizeof(float));
    result->phases          = (float*)malloc(result->num_frequencies * result->output_size * sizeof(float));
    result->phasers         = (float*)malloc(result->num_frequencies * 2 * result->output_size * sizeof(float));

    if (!result->phasers || !result->magnitudes || !result->phases) {
        free_stft(result);
        return;
    }
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

//fftwf : 32bit vs fftw 64bit vs fftwq:  128bit
inline stft_d stft(audio_data *audio,size_t window_size, size_t hop_size,float *window_values) {

    stft_d result      = {0};

    result.sample_rate = audio->sample_rate;

    init_fft_output(&result,window_size,hop_size,audio->num_samples);

    calculate_frequencies(&result,window_size,audio->sample_rate);

    fftwf_complex *in  = (fftwf_complex*)fftwf_malloc(sizeof(fftwf_complex) * window_size);
    fftwf_complex *out = (fftwf_complex*)fftwf_malloc(sizeof(fftwf_complex) * window_size);
    fftwf_plan p       = fftwf_plan_dft_1d((int)window_size, in, out, FFTW_FORWARD, FFTW_ESTIMATE);


    result.output_size = (size_t) result.output_size/audio->channels;


    for (size_t i = 0; i < result.output_size; i++) {
        clock_t start = clock();

        size_t start_idx = i * hop_size * audio->channels; 
        if (start_idx + window_size * audio->channels > audio->num_samples) break;

        for (size_t j = 0; j < window_size; j++) {
            float sum = 0.0f;
            for (int ch = 0; ch < audio->channels; ch++) {
                sum += audio->samples[start_idx + j * audio->channels + ch];
            }

            in[j][0] = (sum / audio->channels) * window_values[j];
            in[j][1] = 0.0;
        }

        fftwf_execute(p);

        clock_t end = clock();

        for (size_t j = 0; j < result.num_frequencies; j++) {
            size_t index = i * (result.num_frequencies * 2) + (j * 2); // Real, imaginary indices

            result.phasers[index]     = out[j][0];
            result.phasers[index + 1] = out[j][1];

            float mag   = sqrt(out[j][0] * out[j][0] + out[j][1] * out[j][1]);
            float phase = atan2(out[j][1], out[j][0]);

            result.magnitudes[i * result.num_frequencies + j] = mag;
            result.phases[i * result.num_frequencies + j]     = phase;
        }

         result.fft_times += ((float)(end - start)) * pow(10, EXP) / CLOCKS_PER_SEC;
    }

    fftwf_destroy_plan(p);
    p = NULL;
    fftwf_free(in);
    in = NULL;

    fftwf_free(out);

    return result;
}

