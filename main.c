#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h> 
#include <omp.h>
#include <string.h>


#define EXP 3  // milliseconds

#include "headers/audio_tools/audio_visualizer.h"

void time_it(struct timeval start, char *str) {
    struct timeval end;
    gettimeofday(&end, NULL);

    long seconds = end.tv_sec - start.tv_sec;
    long microseconds = end.tv_usec - start.tv_usec;
    double elapsed = seconds + microseconds * 1e-6;

    printf("\n%f ms for %s\n", elapsed * 1000.0, str);

    
}

int main(int argc, char *argv[]) {
       if (argc != 10) {
        fprintf(stderr, "Usage: %s <ip_filename> <op_filename> <window_size> <hop_size> <window_type> <number_of_mel_banks> <min_mel> <max_mel> <num_coff>\n", argv[0]);
        return 1;
    }
    

    const char *ip_filename        = argv[1];
    const char *op_filename        = argv[2];
    int window_size                = atoi(argv[3]);
    int hop_size                   = atoi(argv[4]);
    const char *window_type        = argv[5];
    unsigned int num_filters       = (unsigned short int)atoi(argv[6]);
    float min_mel                  = atof(argv[7]);
    float max_mel                  = atof(argv[8]);
    unsigned int num_coff          = (unsigned short int)atoi(argv[9]);
    


    char stft_filename[100], mel_filename[100], mfcc_filename[100]; 

    sprintf(stft_filename, "%s_stft.png", op_filename); 
    sprintf(mel_filename, "%s_mel.png", op_filename);    
    sprintf(mfcc_filename, "%s_mfcc.png", op_filename); 
    
    unsigned char bg_clr[]         = {0,0,0,255};

    audio_data audio               = auto_detect(ip_filename);

    size_t num_frequencies         = (size_t) window_size / 2;

    print_ad(&audio);

    if(audio.samples!=NULL){

        float *window_values           = (float*) malloc(window_size * sizeof(float));    // pre cal to avoid repated calc

        window_function(window_values,window_size,window_type);


        struct timeval start;
        gettimeofday(&start, NULL); 

        stft_d result           = stft(&audio,window_size,hop_size,window_values);
        float *mel_filter_bank  = (float*) calloc((num_frequencies + 1) * (num_filters + 2), sizeof(float)); // type casting to enforce type safety && zero value init ;
        
        cs_from_enum(JET,true);
        
        time_it(start,"STFT");
        gettimeofday(&start, NULL); 
        
        spectrogram(&result, stft_filename,min_mel,max_mel, bg_clr,false,JET);
        time_it(start,"spectrogram");
        gettimeofday(&start, NULL); 
        
        mel_filter(min_mel,max_mel,num_filters,audio.sample_rate,window_size,mel_filter_bank);   // pre cal to avoid repated calc
        float *mel_values = mel_spectrogram(&result,mel_filename,num_filters,mel_filter_bank,bg_clr,true,CIVIDIS);
        time_it(start,"mel spectrogram");
        gettimeofday(&start, NULL); 
        
        mfcc(mel_values,mfcc_filename,result.output_size,num_filters,num_coff,bg_clr,INFERNO);
        time_it(start,"MFCC");

        // #pragma omp parallel sections
        // {
        //     #pragma omp section
        //     {
                 
        //     }
            
        //     #pragma omp section
        //     {
        //         float *mel_values = mel_spectrogram(&result,mel_filename,num_filters,mel_filter_bank,bg_clr,true,1);
        //         mfcc(mel_values,mfcc_filename,result.output_size,num_filters,num_coff,bg_clr,6);
        //     }
        // }
 
        free_stft(&result);

        

    } 
   
    
}
