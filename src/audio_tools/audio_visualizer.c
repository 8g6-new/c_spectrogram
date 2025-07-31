#include "../../headers/audio_tools/audio_visualizer.h"



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


inline float  brachless_db(float mag,bool db){
    return !db*mag + db*log10f(1 + mag*mag);
}

void print_bounds(bounds2d_t *bounds) {
    printf("Time bounds:\n");
    printf("  Start (f): %.2f\n", bounds->time.start_f);
    printf("  End   (f): %.2f\n", bounds->time.end_f);
    printf("  Start (d): %zu\n", bounds->time.start_d);
    printf("  End   (d): %zu\n", bounds->time.end_d);

    printf("Frequency bounds:\n");
    printf("  Start (f): %.2f\n", bounds->freq.start_f);
    printf("  End   (f): %.2f\n", bounds->freq.end_f);
    printf("  Start (d): %zu\n", bounds->freq.start_d);
    printf("  End   (d): %zu\n", bounds->freq.end_d);
}

void set_limits(bounds2d_t *bounds,const size_t max_freq,const size_t max_time){
    bounds->time.end_d  = max_time;
}

void init_bounds(bounds2d_t *bounds,stft_d *result){
    bounds->time.start_d = 0;

    bounds->freq.start_d = hz_to_index(result->num_frequencies,result->sample_rate,bounds->freq.start_f);
    bounds->freq.end_d   = hz_to_index(result->num_frequencies,result->sample_rate,bounds->freq.end_f);
}

void fast_copy(float *data,float *mags,bounds2d_t *bounds,const size_t length){
    const size_t freq_range = bounds->freq.end_d - bounds->freq.start_d;
    const size_t copy_size  = freq_range * sizeof(float);
     
    for (size_t t = bounds->time.start_d; t < bounds->time.end_d; t++) {
        memcpy(data + (t - bounds->time.start_d) * freq_range,mags + t * length + bounds->freq.start_d,copy_size);
    }
}


inline void plot(float *data, bounds2d_t *bounds, plot_t *settings) {
    const size_t w = settings->w;
    const size_t h = settings->h;

    heatmap_t *hm = heatmap_new(w, h);
    if (!hm) {
        fprintf(stderr, "Failed to allocate heatmap.\n");
        return;
    }

    const bool db = settings->db;
    const size_t tstart = bounds->time.start_d;
    const size_t tend = bounds->time.end_d;

    #pragma omp parallel for
    for (size_t t = tstart; t < tend; t++) {
        const size_t offset = (t - tstart) * h;
        for (size_t f = 0; f < h; f++) {
            const float val = brachless_db(data[offset + f], db);
            heatmap_add_weighted_point(hm, t - tstart, h - f - 1, val);
        }
    }

    save_heatmap(&hm, settings->output_file,w, h, settings->bg_color, settings->cs_enum);
}

float *apply_filter_bank(float *cont_mem, size_t num_filters, size_t num_freq, float *mel_filter_bank, bounds2d_t *bounds, plot_t *settings) {
    const bool db       = settings->db;
    const size_t tstart = bounds->time.start_d;
    const size_t tend   = bounds->time.end_d;
    const size_t w      = tend - tstart;
    const size_t h      = bounds->freq.end_d - bounds->freq.start_d;

    float *mel_values   = (float*) malloc(w * h * sizeof(float));
    if (!mel_values) {
        fprintf(stderr, "Failed to allocate mel_values.\n");
        return NULL;
    }

    #pragma omp parallel for
    for (size_t t = tstart; t < tend; t++) {
        size_t offset_in  = (t - tstart) * h;
        size_t offset_out = t * num_filters;

        for (size_t m = 0; m < num_filters; m++) {
            float sum = cblas_sdot(num_freq, &cont_mem[offset_in], 1, &mel_filter_bank[m * num_freq], 1); 
            mel_values[offset_out + m] = brachless_db(sum, db);
        }
    }

  
    return mel_values;
}


float *FCC(float *mel_values, dct_t *dft_coff, bounds2d_t *bounds, plot_t *settings) {
    const bool db        = settings->db;
    const size_t tstart  = bounds->time.start_d;
    const size_t tend    = bounds->time.end_d;
    const size_t w       = tend - tstart;
    const size_t num_f   = dft_coff->num_filters;
    const size_t num_c   = dft_coff->num_coff;

    float *fcc_values = (float*) malloc(w * num_c * sizeof(float));
    if (!fcc_values) {
        fprintf(stderr, "Failed to allocate fcc_values.\n");
        return NULL;
    }

    #pragma omp parallel for
    for (size_t t = tstart; t < tend; t++) {
        size_t offset_in  = (t - tstart) * num_f;
        size_t offset_out = (t - tstart) * num_c;

        for (size_t c = 0; c < num_c; c++) {
            float sum = cblas_sdot(num_f, &dft_coff->coeffs[c * num_f], 1,
                                            &mel_values[offset_in], 1);
            fcc_values[offset_out + c] = brachless_db(sum, db);
        }
    }

    return fcc_values;
}


dct_t gen_cosine_coeffs(const size_t num_filters, const size_t num_coff) {
    dct_t dft_coff = {
        .num_filters = num_filters,
        .num_coff = num_coff,
        .coeffs = malloc(num_filters * num_coff * sizeof(float))
    };

    if (!dft_coff.coeffs) {
        return dft_coff; 
    }

    const float scale = sqrtf(2.0f / (float)num_filters);

    #pragma omp parallel for schedule(static, 256)
    for (size_t n = 0; n < num_coff; n++) {
        for (size_t mel = 0; mel < num_filters; mel++) {
            dft_coff.coeffs[n * num_filters + mel] = scale *
                cosf((float)M_PI / num_filters * (mel + 0.5f) * n);
        }
    }

    return dft_coff;
}
