#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../src/stft/audio_tools.h"

#define EXP 6  // microseconds

float hm_freq = 0.0,std_mean=0.0,std_dev=0.0;

char check(int i,int lim){
    if(i<lim-1) return '\n';
    else return ' ';
}

int main(int argc, char *argv[]) {

    if (argc != 9) {
        fprintf(stderr, "Usage: %s <filename> <window_size> <hop_size> <window_type> <number_of_samples_to_take_average>\n", argv[0]);
        return 1;
    }

    const char *filename            = argv[1];
    int window_size                 = atoi(argv[2]);
    int hop_size                    = atoi(argv[3]);
    const char *window_type         = argv[4];
    unsigned short int num_samples  = (unsigned short int)atoi(argv[5]);
    unsigned short int peak         = (unsigned short int)atoi(argv[6]);
    unsigned short int min          = (unsigned short int)atoi(argv[7]);
    unsigned short int max          = (unsigned short int)atoi(argv[8]);

    double fft_times[window_size];

    if (window_size <= 0 || hop_size <= 0 || num_samples <= 0) {
        fprintf(stderr, "Window size, hop size, and number of samples must be positive integers.\n");
        return 1;
    }


    STFTResult result = stft(filename, window_size, hop_size, window_type,fft_times);
    
    if (result.phasers == NULL || result.frequencies == NULL) {
        fprintf(stderr, "STFT computation failed.\n");
        return 1;
    }
    
    for (int i = 0; i < result.output_size; i++) {
        STFTStats stats = stft_weighted_avg(result.phasers + i * result.num_frequencies * 2, result.num_frequencies, result.frequencies, num_samples,peak,min,max);
        printf("%f,%f,%f,%f,%f,%f,%f,%f,%f%c", 
               /*1, 2, 3, 4, 5, 6, 7, 8, 9*/    stats.mean_freq,
                                                    stats.weighted_mean,
                                                    stats.std_freq_dev,
                                                    fft_times[i],
                                                    stats.time,
                                                    stats.top_freq,
                                                    stats.sdfbpf,
                                                    stats.sdfbmif,
                                                    stats.sdfbmxf,
                                                    check(i, result.output_size)
        );
    }
    

    free(result.phasers);
    free(result.frequencies);

    return 0;
}


