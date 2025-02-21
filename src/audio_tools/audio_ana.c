#include "../../headers/audio_tools/audio_ana.h"


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
    out=NULL;

    free(window_values);
    window_values = NULL;

    return result;
}



int compare_by_magnitude(const void *a, const void *b) {
    float mag_a = ((freq_mag_pair*)a)->magnitude;
    float mag_b = ((freq_mag_pair*)b)->magnitude;
    return (mag_a < mag_b) - (mag_a > mag_b);
}


stft_mean compute_stft_stats(stft_d *result, float f_min, float f_max, bool db) {
    size_t w = result->output_size;
    size_t h = (size_t)(result->num_frequencies * f_max) / (result->sample_rate / 2);
    size_t min_f = (size_t)(result->num_frequencies * f_min) / (result->sample_rate / 2);

    float max = 0.0f;
    float min = 5.0f;

    stft_mean mean = {0};

    size_t size = h - min_f;

    freq_mag_pair *array = calloc(size,sizeof(freq_mag_pair));

    for (size_t t = 0; t < w; t++) {

        for (size_t f = min_f; f < h; f++) {
            float mag = result->magnitudes[t * result->num_frequencies + f];
            mag = db ? mag : log10f(mag * mag + 1.0f);

            if(t==0) array[f - min_f].frequency = result->frequencies[f];  

            array[f - min_f].magnitude += (float) mag/w;
        }
    }

    qsort(array, size, sizeof(freq_mag_pair), compare_by_magnitude);

    for (size_t f = 0; f < size; f++){
        printf("%0.3f %0.3f\n",array[f].frequency,array[f].magnitude/array[0].magnitude);
    }
    

    return mean;
}

