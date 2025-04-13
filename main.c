#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h> 
#include <omp.h>
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
        fprintf(stderr, "Usage: %s <ip_filename> <op_filename> <window_size> <hop_size> <window_type> <number_of_mel_banks> <min_mel> <max_mel> <num_coff> <cs_stft> <cs_mel> <cs_mfcc>\n", argv[0]);
        return 1;
    }
    

    const char *ip_filename        = argv[1];
    const char *op_filename        = argv[2];
    int window_size                = atoi(argv[3]);
    int hop_size                   = atoi(argv[4]);
    const char *window_type        = argv[5];
    unsigned short num_filters     = (unsigned short int)atoi(argv[6]);
    float min_mel                  = atof(argv[7]);
    float max_mel                  = atof(argv[8]);
    unsigned short num_coff        = (unsigned short int)atoi(argv[9]);
    unsigned short cs_stft         = (unsigned short int)atoi(argv[10]);
    unsigned short cs_mel          = (unsigned short int)atoi(argv[11]);
    unsigned short cs_mfcc         = (unsigned short int)atoi(argv[12]);
    const char         *cache_fol  = argv[13];

    char stft_filename[100], mel_filename[100], mfcc_filename[100]; 

    sprintf(stft_filename, "%s_stft.png", op_filename); 
    sprintf(mel_filename, "%s_mel.png", op_filename);    
    sprintf(mfcc_filename, "%s_mfcc.png", op_filename); 
    
    unsigned char bg_clr[]         = {0,0,0,255};
    
    START_TIMING();
    audio_data audio               = auto_detect(ip_filename);
    END_TIMING("auto_det");

    size_t num_frequencies         = (size_t) window_size / 2;

    print_ad(&audio);

    if(audio.samples!=NULL){


        START_TIMING();
        float *window_values           = (float*) malloc(window_size * sizeof(float));    // pre cal to avoid repated calc
        window_function(window_values,window_size,window_type);
        float *mel_filter_bank  = (float*) calloc((num_frequencies + 1) * (num_filters + 2), sizeof(float)); // type casting to enforce type safety && zero value init ;
        mel_filter(min_mel,max_mel,num_filters,audio.sample_rate,window_size,mel_filter_bank);   // pre cal to avoid repated calc
        END_TIMING("init");


        START_TIMING();
        fft_d fft_plan = init_fftw_plan(window_size,cache_fol);
        END_TIMING("fetch fft cache 1");


        cs_from_enum(cs_stft,true);

        START_TIMING();
        stft_d result           = stft(&audio,window_size,hop_size,window_values,&fft_plan);
        END_TIMING("stft");


        START_TIMING();
        spectrogram(&result, stft_filename,min_mel,max_mel, bg_clr,false,cs_stft);
        END_TIMING("stft_plot");

        
        START_TIMING();
        float *mel_values = mel_spectrogram(&result,mel_filename,num_filters,mel_filter_bank,bg_clr,true,cs_mel);
        END_TIMING("mel");

        START_TIMING();
        mfcc(mel_values,mfcc_filename,result.output_size,num_filters,num_coff,bg_clr,cs_mfcc);
        END_TIMING("mfcc");

        print_bench_ranked();
 
        free_stft(&result);

    } 
   
    
}
