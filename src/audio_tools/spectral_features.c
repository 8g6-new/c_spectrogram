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

#define CALCULATE_MAGNITUDE 1
#define CALCULATE_PHASE 1

inline size_t hz_to_index(size_t num_freq, size_t sample_rate, float f) {
    return (size_t)((num_freq * f) / (sample_rate / 2.0f));
}

inline size_t align32(size_t size) {
    return (size + 31) & ~((size_t)31);
}

static inline void compute_magnitudes_and_phases_scalar(fftwf_complex *out, float *magnitudes, float *phases, const size_t offset, const size_t count) {
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


static inline __m128 atan2_ps_sse(__m128 y, __m128 x) {
    const __m128 zero = _mm_set1_ps(0.0f);
    const __m128 one = _mm_set1_ps(1.0f);
    const __m128 pi = _mm_set1_ps(3.14159265358979323846f);
    const __m128 pi_2 = _mm_set1_ps(1.57079632679489661923f);
    const __m128 pi_4 = _mm_set1_ps(0.785398163397448309616f);
    const __m128 c = _mm_set1_ps(0.273f);

    // Prevent divide-by-zero
    const __m128 safe_x = _mm_blendv_ps(x, _mm_set1_ps(1e-10f), _mm_cmpeq_ps(x, zero));
    const __m128 z = _mm_div_ps(y, safe_x);
    const __m128 abs_z = _mm_andnot_ps(_mm_set1_ps(-0.0f), z); // clear sign bit

    // atan(z) ≈ z * pi/4 + c * z * (|z| - 1)
    const __m128 term = _mm_mul_ps(c, _mm_mul_ps(_mm_sub_ps(abs_z, one), z));
    const __m128 atan = _mm_fmadd_ps(z, pi_4, term);

    // Adjust results based on quadrant
    __m128 result = atan;
    const __m128 x_neg = _mm_cmplt_ps(x, zero);
    const __m128 y_neg = _mm_cmplt_ps(y, zero);
    const __m128 y_pos = _mm_cmpgt_ps(y, zero);
    const __m128 x_zero = _mm_cmpeq_ps(x, zero);

    result = _mm_blendv_ps(result, _mm_add_ps(atan, pi), _mm_andnot_ps(y_neg, x_neg));      // x < 0, y ≥ 0
    result = _mm_blendv_ps(result, _mm_sub_ps(atan, pi), _mm_and_ps(x_neg, y_neg));          // x < 0, y < 0
    result = _mm_blendv_ps(result, pi_2, _mm_and_ps(x_zero, y_pos));                            // x == 0, y > 0
    result = _mm_blendv_ps(result, _mm_sub_ps(zero, pi_2), _mm_and_ps(x_zero, y_neg));       // x == 0, y < 0

    return result;
}

static inline void compute_magnitudes_and_phases_sse(fftwf_complex *out, float *magnitudes, float *phases, const size_t offset, const size_t num_frequencies) {
    const size_t simd_limit = num_frequencies & ~3UL; // multiple of 4
    float *real_ptr = &out[0][0];
    float *imag_ptr = &out[0][1];
    const size_t stride = 2;

    for (size_t j = 0; j < simd_limit; j += 4) {
        __m128 real_vals = _mm_loadu_ps(real_ptr + j * stride);
        __m128 imag_vals = _mm_loadu_ps(imag_ptr + j * stride);

        #ifdef CALCULATE_MAGNITUDE
                __m128 imag_squared = _mm_mul_ps(imag_vals, imag_vals);
                __m128 sum = _mm_add_ps(_mm_mul_ps(real_vals, real_vals), imag_squared);
                __m128 magnitudes_vals = _mm_sqrt_ps(sum);
                _mm_storeu_ps(&magnitudes[offset + j], magnitudes_vals);
        #endif

        #ifdef CALCULATE_PHASE
                __m128 phase_vals = atan2_ps_sse(imag_vals, real_vals);
                _mm_storeu_ps(&phases[offset + j], phase_vals);
        #endif
    }

    if (simd_limit < num_frequencies) {
        compute_magnitudes_and_phases_scalar(out + simd_limit, magnitudes, phases, offset + simd_limit, num_frequencies - simd_limit);
    }
}

static inline __m256 atan2_ps_avx(__m256 y, __m256 x) {
    const __m256 zero = _mm256_set1_ps(0.0f);
    const __m256 one = _mm256_set1_ps(1.0f);
    const __m256 pi = _mm256_set1_ps(3.14159265358979323846f);
    const __m256 pi_2 = _mm256_set1_ps(1.57079632679489661923f);
    const __m256 pi_4 = _mm256_set1_ps(0.785398163397448309616f);
    const __m256 c = _mm256_set1_ps(0.273f);

    // Prevent divide-by-zero
    const __m256 safe_x = _mm256_blendv_ps(x, _mm256_set1_ps(1e-10f), _mm256_cmp_ps(x, zero, _CMP_EQ_OS));
    const __m256 z = _mm256_div_ps(y, safe_x);
    const __m256 abs_z = _mm256_andnot_ps(_mm256_set1_ps(-0.0f), z); // clear sign bit

    // atan(z) ≈ z * pi/4 + c * z * (|z| - 1)
    const __m256 term = _mm256_mul_ps(c, _mm256_mul_ps(_mm256_sub_ps(abs_z, one), z));
    const __m256 atan = _mm256_fmadd_ps(z, pi_4, term);

    // Adjust results based on quadrant
    __m256 result = atan;
    const __m256 x_neg = _mm256_cmp_ps(x, zero, _CMP_LT_OS);
    const __m256 y_neg = _mm256_cmp_ps(y, zero, _CMP_LT_OS);
    const __m256 y_pos = _mm256_cmp_ps(y, zero, _CMP_GT_OS);
    const __m256 x_zero = _mm256_cmp_ps(x, zero, _CMP_EQ_OS);

    result = _mm256_blendv_ps(result, _mm256_add_ps(atan, pi), _mm256_andnot_ps(y_neg, x_neg));      // x < 0, y ≥ 0
    result = _mm256_blendv_ps(result, _mm256_sub_ps(atan, pi), _mm256_and_ps(x_neg, y_neg));          // x < 0, y < 0
    result = _mm256_blendv_ps(result, pi_2, _mm256_and_ps(x_zero, y_pos));                            // x == 0, y > 0
    result = _mm256_blendv_ps(result, _mm256_sub_ps(zero, pi_2), _mm256_and_ps(x_zero, y_neg));       // x == 0, y < 0

    return result;
}

static inline void compute_magnitudes_and_phases_avx(fftwf_complex *out, float *magnitudes, float *phases, const size_t offset, const size_t num_frequencies) {
    const size_t simd_limit = num_frequencies & ~7UL; // multiple of 8
    float *real_ptr = &out[0][0];
    float *imag_ptr = &out[0][1];
    const size_t stride = 2;

    for (size_t j = 0; j < simd_limit; j += 8) {
        __m256 real_vals = _mm256_loadu_ps(real_ptr + j * stride);
        __m256 imag_vals = _mm256_loadu_ps(imag_ptr + j * stride);

        #ifdef CALCULATE_MAGNITUDE
        // FMA: real*real + imag*imag
        __m256 imag_squared = _mm256_mul_ps(imag_vals, imag_vals);
        __m256 sum = _mm256_fmadd_ps(real_vals, real_vals, imag_squared);
        __m256 magnitudes_vals = _mm256_sqrt_ps(sum);
        _mm256_storeu_ps(&magnitudes[offset + j], magnitudes_vals);
        #endif

        #ifdef CALCULATE_PHASE
        __m256 phase_vals = atan2_ps_avx(imag_vals, real_vals); // your own AVX atan2
        _mm256_storeu_ps(&phases[offset + j], phase_vals);
        #endif
    }

    if (simd_limit < num_frequencies) {
        compute_magnitudes_and_phases_scalar(out + simd_limit, magnitudes, phases, offset + simd_limit, num_frequencies - simd_limit);
    }
}




static inline void compute_magnitudes_and_phases_sse2(fftwf_complex *out, float *magnitudes, float *phases, const size_t offset, const size_t num_frequencies) {
    const size_t simd_limit = num_frequencies & ~3UL; // Round down to multiple of 4
    float *real_ptr = &out[0][0];
    float *imag_ptr = &out[0][1];
    const size_t stride = 2;

    for (size_t j = 0; j < simd_limit; j += 4) {
        __m128 real_vals = _mm_loadu_ps(real_ptr + j * stride);
        __m128 imag_vals = _mm_loadu_ps(imag_ptr + j * stride);

        #ifdef CALCULATE_MAGNITUDE
        __m128 real_squared = _mm_mul_ps(real_vals, real_vals);
        __m128 imag_squared = _mm_mul_ps(imag_vals, imag_vals);
        __m128 sum = _mm_add_ps(real_squared, imag_squared);
        __m128 magnitudes_vals = _mm_sqrt_ps(sum);
        _mm_storeu_ps(&magnitudes[offset + j], magnitudes_vals);
        #endif

        #ifdef CALCULATE_PHASE
        __m128 phase_vals = atan2_ps_sse(imag_vals, real_vals);
        _mm_storeu_ps(&phases[offset + j], phase_vals);
        #endif
    }

    if (simd_limit < num_frequencies) {
        compute_magnitudes_and_phases_scalar(out + simd_limit, magnitudes, phases, offset + simd_limit, num_frequencies - simd_limit);
    }
}

// static inline __m512 atan2_ps_avx512(__m512 y, __m512 x) {
//     const __m512 zero = _mm512_set1_ps(0.0f);
//     const __m512 one = _mm512_set1_ps(1.0f);
//     const __m512 pi = _mm512_set1_ps(3.14159265358979323846f);
//     const __m512 pi_2 = _mm512_set1_ps(1.57079632679489661923f);
//     const __m512 pi_4 = _mm512_set1_ps(0.785398163397448309616f);
//     const __m512 c = _mm512_set1_ps(0.273f);

//     // Prevent divide-by-zero
//     const __m512 safe_x = _mm512_mask_blend_ps(_mm512_cmpeq_ps_mask(x, zero), x, _mm512_set1_ps(1e-10f));
//     const __m512 z = _mm512_div_ps(y, safe_x);
//     const __m512 abs_z = _mm512_andnot_ps(_mm512_set1_ps(-0.0f), z); // clear sign bit

//     // atan(z) ≈ z * pi/4 + c * z * (|z| - 1)
//     const __m512 term = _mm512_mul_ps(c, _mm512_mul_ps(_mm512_sub_ps(abs_z, one), z));
//     const __m512 atan = _mm512_fmadd_ps(z, pi_4, term);

//     // Adjust results based on quadrant
//     __m512 result = atan;
//     const __mmask16 x_neg_mask = _mm512_cmplt_ps_mask(x, zero);
//     const __mmask16 y_neg_mask = _mm512_cmplt_ps_mask(y, zero);
//     const __mmask16 y_pos_mask = _mm512_cmpgt_ps_mask(y, zero);
//     const __mmask16 x_zero_mask = _mm512_cmpeq_ps_mask(x, zero);

//     result = _mm512_mask_blend_ps(_mm512_kandn(y_neg_mask, x_neg_mask), result, _mm512_add_ps(atan, pi));     // x < 0, y ≥ 0
//     result = _mm512_mask_blend_ps(_mm512_kand(x_neg_mask, y_neg_mask), result, _mm512_sub_ps(atan, pi));       // x < 0, y < 0
//     result = _mm512_mask_blend_ps(_mm512_kand(x_zero_mask, y_pos_mask), result, pi_2);                         // x == 0, y > 0
//     result = _mm512_mask_blend_ps(_mm512_kand(x_zero_mask, y_neg_mask), result, _mm512_sub_ps(zero, pi_2));    // x == 0, y < 0

//     return result;
// }

// static inline  void compute_magnitudes_and_phases_avx512(fftwf_complex *out, float *magnitudes, float *phases, const size_t offset, const size_t num_frequencies) {
//     const size_t simd_limit = num_frequencies & ~15UL; // multiple of 16
//     float *real_ptr = &out[0][0];
//     float *imag_ptr = &out[0][1];
//     const size_t stride = 2;

//     for (size_t j = 0; j < simd_limit; j += 16) {
//         __m512 real_vals = _mm512_loadu_ps(real_ptr + j * stride);
//         __m512 imag_vals = _mm512_loadu_ps(imag_ptr + j * stride);

//         #ifdef CALCULATE_MAGNITUDE
//                 __m512 imag_squared = _mm512_mul_ps(imag_vals, imag_vals);
//                 __m512 sum = _mm512_fmadd_ps(real_vals, real_vals, imag_squared);
//                 __m512 magnitudes_vals = _mm512_sqrt_ps(sum);
//                 _mm512_storeu_ps(&magnitudes[offset + j], magnitudes_vals);
//         #endif

//         #ifdef CALCULATE_PHASE
//                 __m512 phase_vals = atan2_ps_avx512(imag_vals, real_vals);
//                 _mm512_storeu_ps(&phases[offset + j], phase_vals);
//         #endif
//     }

//     if (simd_limit < num_frequencies) {
//         compute_magnitudes_and_phases_scalar(out + simd_limit, magnitudes, phases, offset + simd_limit, num_frequencies - simd_limit);
//     }
// }

// void compute_magnitudes_and_phases_simd(fftwf_complex *out, float *magnitudes, float *phases, const size_t offset, const size_t num_frequencies) {
//     #if defined(__AVX512F__)
//         compute_magnitudes_and_phases_avx512(out, magnitudes, phases, offset, num_frequencies);
//     #elif defined(__AVX__)
//         compute_magnitudes_and_phases_avx(out, magnitudes, phases, offset, num_frequencies);
//     #elif defined(__SSE2__)
//         compute_magnitudes_and_phases_sse2(out, magnitudes, phases, offset, num_frequencies);
//     #else
//         compute_magnitudes_and_phases_scalar(out, magnitudes, phases, offset, num_frequencies);
//     #endif
// }





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

    if(hop_size<=0)
      perror("Hop size should be >= 0");

    if(window_size<=0)
      perror("Hop size should be >= 0");

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

    result->magnitudes      = (float*)aligned_alloc(32, align32(magnitudes_size));
    result->phases          = (float*)aligned_alloc(32, align32(phases_size));
    result->phasers         = (float*)aligned_alloc(32, align32(phasers_size));


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

float safe_div(float a,float b){
    return b==0 ? 0.0f : a/b;
}

melbank_t mel_filter(float min_f, float max_f, size_t n_filters, float sr, size_t fft_size, float *filter) {

    melbank_t non_zero      = { .freq_indexs = NULL, .weights = NULL, .size = 0 ,.num_filters=n_filters };
   
    const size_t avg_length = n_filters * (fft_size * (fft_size*5)/512);
    non_zero.freq_indexs    = malloc(avg_length * sizeof(size_t));
    non_zero.weights        = malloc(avg_length * sizeof(float));
    
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
                    (num_samples - window_size * channels) / hop_size : 0;
    max_i = (max_i < output_size) ? max_i : output_size;

    result.output_size = max_i;

    // Allocate one FFT context
    fftwf_complex* in  = (fftwf_complex*)aligned_alloc(32, align32(window_size * sizeof(fftwf_complex)));
    fftwf_complex* out = (fftwf_complex*)aligned_alloc(32, align32(window_size * sizeof(fftwf_complex)));

    fftwf_plan plan = fftwf_plan_dft_1d(
        (int)window_size, 
        in, 
        out, 
        FFTW_FORWARD, 
        FFTW_WISDOM_ONLY
    );

    if (!plan) {
        plan = fftwf_plan_dft_1d(
            (int)window_size, 
            in, 
            out, 
            FFTW_FORWARD, 
            FFTW_ESTIMATE
        );
    }

    for (size_t i = 0; i < max_i; i++) {
        const size_t start_idx = i * hop_size;

        for (size_t j = 0; j < window_size; j++) {
            float sum = 0.0f;
            for (size_t ch = 0; ch < channels; ch++) {
                sum += audio->samples[start_idx + (j * channels) + ch];
            }
            in[j][0] = (sum * window_values[j]) / channels;
            in[j][1] = 0.0f;
        }

        fftwf_execute(plan);

        const size_t offset = i * num_frequencies;
        compute_magnitudes_and_phases_scalar(out, result.magnitudes, result.phases, offset, num_frequencies);
    }

    fftwf_destroy_plan(plan);
    free(in);
    free(out);

    return result;
}
