#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <omp.h>

#include "src/stft/audio_tools.h"

#define EXP 6  // microseconds


int main(int argc, char *argv[]) {
    if (argc != 8) {
        fprintf(stderr, "Usage: %s <ip_filename> <op_filename> <window_size> <hop_size> <window_type> <number_of_mel_banks> <num_coff>\n", argv[0]);
        return 1;
    }
    
    clock_t start = clock();

    const char *ip_filename        = argv[1];
    const char *op_filename        = argv[2];
    int window_size                = atoi(argv[3]);
    int hop_size                   = atoi(argv[4]);
    const char *window_type        = argv[5];
    unsigned int num_filters       = (unsigned short int)atoi(argv[6]);
    unsigned int num_coff          = (unsigned short int)atoi(argv[7]);


    char stft_filename[100], mel_filename[100], mfcc_filename[100]; 

    sprintf(stft_filename, "%s_stft.png", op_filename); 
    sprintf(mel_filename, "%s_mel.png", op_filename);    
    sprintf(mfcc_filename, "%s_mfcc.png", op_filename); 
    
    unsigned char bg_clr[]         = {0,0,0,255};

    audio_data audio               = read_file(ip_filename);

    size_t num_frequencies         = (size_t) window_size / 2; 

    float *frequencies             = (float*) malloc(num_frequencies * sizeof(float));
    float *window_values           = (float*) malloc(window_size * sizeof(float));

    calculate_frequencies(frequencies, window_size,audio.sample_rate);
    window_function(window_values,window_size,window_type);

    stft_d result = stft(
            audio.samples,
            window_size,hop_size,audio.num_samples,window_values
    );
    
    float *mel_values              = (float*) malloc(num_filters * result.output_size * sizeof(float));
 
    size_t w                       = result.output_size;


    float *mel_filter_bank    = (float*) calloc((num_frequencies + 1) * (num_filters + 2), sizeof(float)); // type casting to enforce type safety && zero value init ;
    mel_filter(0.0,audio.sample_rate/2,num_filters,audio.sample_rate,window_size,mel_filter_bank);
    
    #pragma omp parallel sections
    {
        #pragma omp section
        {
            spectrogram(&result, stft_filename, bg_clr, true);
        }
        
        #pragma omp section
        {
            mel_spectrogram(&result,mel_filename, bg_clr, true, num_filters, mel_filter_bank, mel_values);
            MFCC(mel_values,mfcc_filename, w, bg_clr, true, num_filters, num_coff);
        }
    }
 
    
    free_stft(&result);

    clock_t end = clock();
    printf("%f ms\n",(float)(end - start) * pow(10,3) / CLOCKS_PER_SEC);
}
