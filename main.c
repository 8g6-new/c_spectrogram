#include <stdio.h>
#include <stdlib.h>
#include <string.h>


#include "headers/audio_tools/audio_visualizer.h"
#include "headers/utils/bench.h"



int main(int argc, char *argv[]) {
       if (argc != 14) {
        fprintf(stderr, "Usage: %s <ip_filename> <op_filename> <window_size> <hop_size> <window_type> <number_of_mel_banks> <min_mel> <max_mel> <num_coff> <cs_stft> <cs_mel> <cs_mfcc> <cache_fol>\n", argv[0]);
        return 1;
    }
    

    const char *ip_filename        = argv[1];
    const char *op_filename        = argv[2];
    int window_size                = atoi(argv[3]);
    int hop_size                   = atoi(argv[4]);
    const char *window_type        = argv[5];
    const size_t   num_filters     = (unsigned short int)atoi(argv[6]);
    float min_mel                  = atof(argv[7]);
    float max_mel                  = atof(argv[8]);
    unsigned short num_coff        = (unsigned short int)atoi(argv[9]);
    unsigned short cs_stft         = (unsigned short int)atoi(argv[10]);
    unsigned short cs_mel          = (unsigned short int)atoi(argv[11]);
    unsigned short cs_mfcc         = (unsigned short int)atoi(argv[12]);
    const char         *cache_fol  = argv[13];

    char *stft_str = cs_from_enum(cs_stft, false);
    char *mel_str  = cs_from_enum(cs_mel, false);
    char *mfcc_str = cs_from_enum(cs_mfcc, false);



    printf("Input Filename       : %s\n", ip_filename);
    printf("Output Filename      : %s\n", op_filename);
    printf("Window Size          : %d\n", window_size);
    printf("Hop Size             : %d\n", hop_size);
    printf("Window Type          : %s\n", window_type);
    printf("Number of Filters    : %zu\n", num_filters);
    printf("Min Mel Frequency    : %.2f\n", min_mel);
    printf("Max Mel Frequency    : %.2f\n", max_mel);
    printf("Number of Coeffs     : %hu\n", num_coff);
    printf("STFT Color Scheme    : %s\n", stft_str);
    printf("Mel Color Scheme     : %s\n", mel_str);
    printf("MFCC Color Scheme    : %s\n", mfcc_str);
    printf("Cache Folder         : %s\n", cache_fol);

    free(stft_str);
    free(mel_str);
    free(mfcc_str);
    
    audio_data audio               = auto_detect(ip_filename);

    size_t num_frequencies         = (size_t) window_size / 2;

    print_ad(&audio);

    if (audio.samples != NULL) {

        // ---- Precompute: Window ----
        START_TIMING();
        float *window_values = (float*) malloc(window_size * sizeof(float));
        window_function(window_values, window_size, window_type);
        END_TIMING("pre:win");

        // ---- Precompute: Mel Filterbank ----
        START_TIMING();
        float *filterbank = (float*) calloc((num_frequencies + 1) * (num_filters + 2), sizeof(float));
        filter_bank_t bank = gen_filterbank(F_MEL, min_mel, max_mel, num_filters,
                                            audio.sample_rate, window_size, filterbank);
        END_TIMING("pre:mel");

        // ---- Precompute: DCT Coefficients ----
        START_TIMING();
        dct_t dft_coff = gen_cosine_coeffs(num_filters, num_coff);
        END_TIMING("pre:dct");

        // ---- Precompute: FFT Plan ----
        START_TIMING();
        fft_d fft_plan = init_fftw_plan(window_size, cache_fol);
        END_TIMING("pre:fft");

        // ---- Plot Settings ----
        cs_from_enum(cs_stft, true);
        plot_t settings = {
            .cs_enum = cs_stft,
            .db = true
        };
        settings.bg_color[0] = 0;
        settings.bg_color[1] = 0;
        settings.bg_color[2] = 0;
        settings.bg_color[3] = 255;

        // ---- STFT ----
        START_TIMING();
        stft_d result = stft(&audio, window_size, hop_size, window_values, &fft_plan);
        END_TIMING("stft");

        printf("\nstft ended\n");

        // ---- Bounds ----
        bounds2d_t bounds = {0};
        bounds.freq.start_f = min_mel;
        bounds.freq.end_f   = max_mel;

        init_bounds(&bounds, &result);
        set_limits(&bounds, result.num_frequencies, result.output_size);

        print_bounds(&bounds);

        const size_t t_len = bounds.time.end_d - bounds.time.start_d;
        const size_t f_len = bounds.freq.end_d - bounds.freq.start_d;

        float *cont_mem = malloc(t_len * f_len * sizeof(float));

        // ---- Copy Magnitudes ----
        START_TIMING();
        fast_copy(cont_mem, result.magnitudes, &bounds, result.num_frequencies);
        END_TIMING("copy");

        printf("\ncopy ended\n");

        // ---- Plot STFT ----
        sprintf(settings.output_file, "%s_stft.png", op_filename);
        settings.w = t_len;
        settings.h = f_len;
        START_TIMING();
        plot(cont_mem, &bounds, &settings);
        END_TIMING("plt:stft");

        // ---- Compute Mel ----
        START_TIMING();
        float *mel_values = apply_filter_bank(cont_mem, num_filters, result.num_frequencies,
                                            filterbank, &bounds, &settings);
        END_TIMING("mel");

        // ---- Plot Mel ----
        sprintf(settings.output_file, "%s_mel.png", op_filename);
        settings.h = num_filters;
        settings.w = t_len;
        START_TIMING();
        plot(mel_values, &bounds, &settings);
        END_TIMING("plt:mel");


        // ---- Compute MFCC ----
        START_TIMING();
        float *fcc_values = FCC(mel_values, &dft_coff, &bounds, &settings);
        END_TIMING("mfcc");

        // ---- Plot MFCC ----
        sprintf(settings.output_file, "%s_mfcc.png", op_filename);
        settings.h = dft_coff.num_coff;
        settings.w = t_len;
        START_TIMING();
        plot(fcc_values, &bounds, &settings);
        END_TIMING("plt:mfcc");

        // ---- Cleanup ----
        free_stft(&result);
        free_audio(&audio);
        free(window_values);
        free(mel_values);
        free(fcc_values);
        free(filterbank);
        free(cont_mem);
        free_fft_plan(&fft_plan);
        free(dft_coff.coeffs);
        free(bank.freq_indexs);
        free(bank.weights);

        print_bench_ranked();
    }

    
}
