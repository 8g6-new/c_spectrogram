#include <stdio.h>
#include <stdlib.h>
#include <string.h>


#define EXP 3  // milliseconds

#include "headers/audio_tools/audio_visualizer.h"
#include "headers/utils/bench.h"

void time_it(struct timeval start, char *str) {
    struct timeval end;
    gettimeofday(&end, NULL);

    long seconds = end.tv_sec - start.tv_sec;
    long microseconds = end.tv_usec - start.tv_usec;
    double elapsed = seconds + microseconds * 1e-6;

    printf("\n%f ms for %s\n", elapsed * 1000.0, str);

    
}

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


    printf("Input Filename       : %s\n", ip_filename);
    printf("Output Filename      : %s\n", op_filename);
    printf("Window Size          : %d\n", window_size);
    printf("Hop Size             : %d\n", hop_size);
    printf("Window Type          : %s\n", window_type);
    printf("Number of Filters    : %zu\n", num_filters);
    printf("Min Mel Frequency    : %.2f\n", min_mel);
    printf("Max Mel Frequency    : %.2f\n", max_mel);
    printf("Number of Coeffs     : %hu\n", num_coff);
    printf("Cache STFT Channels  : %hu\n", cs_stft);
    printf("Cache Mel Channels   : %hu\n", cs_mel);
    printf("Cache MFCC Channels  : %hu\n", cs_mfcc);
    printf("Cache Folder         : %s\n", cache_fol);


    char stft_filename[100], mel_filename[100], mfcc_filename[100]; 

    audio_data audio               = auto_detect(ip_filename);

    size_t num_frequencies         = (size_t) window_size / 2;

    print_ad(&audio);

    if(audio.samples!=NULL){



        START_TIMING();
        float *window_values           = (float*) malloc(window_size * sizeof(float));    // pre cal to avoid repated calc
        window_function(window_values,window_size,window_type);
        END_TIMING("stft window");
        START_TIMING();
        START_TIMING();
        float *mel_filter_bank  = (float*) calloc((num_frequencies + 1) * (num_filters + 2), sizeof(float)); // type casting to enforce type safety && zero value init ;
        melbank_t non_zero      = mel_filter(min_mel,max_mel,num_filters,audio.sample_rate,window_size,mel_filter_bank);   // pre cal to avoid repated calc
        END_TIMING("mel filter bank");
        START_TIMING();
        mffc_t dft_coff         = precompute_cosine_coeffs(num_filters,num_coff);
        END_TIMING("dtft coff");

    

        START_TIMING();
        fft_d fft_plan = init_fftw_plan(window_size,cache_fol);
        END_TIMING("fetch fft cache 1");

        cs_from_enum(cs_stft,true);

        plot_t settings={
            .cs_enum = cs_stft,
            .db = true
        };

        settings.bg_color[0] = 0;
        settings.bg_color[1] = 0;
        settings.bg_color[2] = 0;
        settings.bg_color[3] = 255;



        START_TIMING();
        stft_d result           = stft(&audio,window_size,hop_size,window_values,&fft_plan);
        END_TIMING("stft");

        printf("\nstft ended\n");

        bounds2d_t bounds   = {0};
        
        bounds.freq.start_f = min_mel;
        bounds.freq.end_f   = max_mel;
        
        init_bounds(&bounds,&result);
        set_limits(&bounds,result.num_frequencies,result.output_size);
       
       
        print_bounds(&bounds);
        float *contious_mem     = malloc((bounds.freq.end_d - bounds.freq.start_d) * (bounds.time.end_d - bounds.time.start_d) * sizeof(float));
        

        START_TIMING();
        fast_copy(contious_mem,result.magnitudes,&bounds,result.num_frequencies);
        END_TIMING("copy");
        printf("\ncopy ended\n");

        free_stft(&result);
        free_audio(&audio);
        
        sprintf(settings.output_file, "%s_stft.png", op_filename); 
        
        START_TIMING();
        spectrogram(contious_mem,&bounds,&settings);
        END_TIMING("stft_plot");

       

        sprintf(settings.output_file, "%s_mel.png", op_filename); 
        START_TIMING();
        float *mel_values = mel_spectrogram(contious_mem,num_filters,result.num_frequencies,mel_filter_bank,&bounds,&settings);
        END_TIMING("mel");
        
        sprintf(settings.output_file, "%s_mfcc.png", op_filename); 
        START_TIMING();
        mfcc(mel_values,&dft_coff,&bounds,&settings);
        END_TIMING("mfcc");

        print_bench_ranked();
 
       
        free(window_values);
        free(mel_values);
        free(mel_filter_bank);
        free(contious_mem);
        free_fft_plan(&fft_plan);
        free(dft_coff.coeffs);
        free(non_zero.freq_indexs);
        free(non_zero.weights);

    } 
   
    
}
