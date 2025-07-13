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

void fast_copy(float *contious_mem,float *mags,bounds2d_t *bounds,const size_t length){
    const size_t freq_range = bounds->freq.end_d - bounds->freq.start_d;
    const size_t copy_size  = freq_range * sizeof(float);
     
    for (size_t t = bounds->time.start_d; t < bounds->time.end_d; t++) {
        memcpy(contious_mem + (t - bounds->time.start_d) * freq_range,mags + t * length + bounds->freq.start_d,copy_size);
    }
}


inline void spectrogram(float *contious_mem, bounds2d_t *bounds, plot_t *settings) {

    
    const size_t w       = bounds->time.end_d - bounds->time.start_d;
    const size_t h       = bounds->freq.end_d - bounds->freq.start_d;

    heatmap_t *hm        = heatmap_new(w,h);
    
    if (!hm) {
        fprintf(stderr, "Failed to allocate heatmap.\n");
    }

    const bool db        = settings->db;
    const size_t tstart  = bounds->time.start_d;
    const size_t tend    = bounds->time.end_d;

     
    #pragma omp parallel for 
    for (size_t t = tstart; t < tend; t++) {
        const size_t offset = (t - tstart) * h;
        for (size_t f = 0; f < h; f++) {
            const float val = brachless_db(contious_mem[offset + f], db);
            heatmap_add_weighted_point(hm, t - tstart, h-f-1, val);
        }
    }

    save_heatmap(&hm, settings->output_file, w, h, settings->bg_color, settings->cs_enum);
}

float *mel_spectrogram(float *contious_mem,const size_t num_filters,const size_t num_freq,float *mel_filter_bank, bounds2d_t *bounds, plot_t *settings){


    omp_set_num_threads(omp_get_max_threads());

    const size_t w       = bounds->time.end_d - bounds->time.start_d;
    const size_t h       = bounds->freq.end_d - bounds->freq.start_d;
    const bool db        = settings->db;
    const size_t tstart  = bounds->time.start_d;
    const size_t tend    = bounds->time.end_d;

    heatmap_t *hm        = heatmap_new(w,num_filters);
    
    if (!hm) {
        fprintf(stderr, "Failed to allocate heatmap.\n");
    }

    float *mel_values    = (float*) malloc(num_filters * w * sizeof(float));

    #pragma omp parallel for she
    for (size_t t = tstart; t < tend; t++) {
        const size_t offset1 = (t - tstart) * h;
        const size_t offset3 = t * num_filters;
        float sum;
        for (size_t mel = 0; mel < num_filters; mel++) {
            const size_t offset2      = mel * num_freq;
            sum                       = cblas_sdot(h, &contious_mem[offset1], 1, &mel_filter_bank[offset2], 1);
            sum                       = brachless_db(sum,db);
            mel_values[offset3 + mel] = sum;
            heatmap_add_weighted_point(hm, t - tstart, num_filters - mel, sum);
          
        }
    }

    
    save_heatmap(&hm, settings->output_file,w,num_filters,settings->bg_color, settings->cs_enum);
   
    return mel_values;
}



mffc_t precompute_cosine_coeffs(const size_t num_filters, const size_t num_coff) {
    mffc_t dft_coff = {
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

inline void mfcc(float *mel_values, mffc_t *dft_coff, bounds2d_t *bounds, plot_t *settings) {
    const bool db            = settings->db;
    const size_t tstart      = bounds->time.start_d;
    const size_t tend        = bounds->time.end_d;
    const size_t w           = tend - tstart;
    const size_t num_filters = dft_coff->num_filters;
    const size_t num_coff    = dft_coff->num_coff;

    heatmap_t *hm            = heatmap_new( w, num_coff);

    if (!hm) {
        fprintf(stderr, "Failed to allocate heatmap.\n");
    }


    #pragma omp parallel for 
    for (size_t t = tstart; t < tend; t++) {
        for (size_t n = 0; n < num_coff; n++) {
            float sum = cblas_sdot(num_filters,&dft_coff->coeffs[n * num_filters], 1,&mel_values[t * num_filters], 1);
            heatmap_add_weighted_point(hm, t - tstart,num_coff - n - 1,brachless_db(sum, db));
        }
    }

    save_heatmap(&hm, settings->output_file, w, num_coff, settings->bg_color, settings->cs_enum);
}
