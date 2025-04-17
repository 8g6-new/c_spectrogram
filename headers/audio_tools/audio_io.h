#ifndef AUDIO_IO_H
#define AUDIO_IO_H

    #define EXP 3  // milliseconds

    #include <stdio.h>       
    #include <stdlib.h>    
    #include <stdint.h>      
    #include <sys/stat.h>    
    #include <sndfile.h>
    #include <time.h>
    #include <math.h>
    
    #include "../utils/ftype_detect.h"
    // #include "../../headers/ml/attention.h"
    


    #define MINIMP3_ONLY_MP3
    #define MINIMP3_USE_SIMD
    #define MINIMP3_IMPLEMENTATION


    #ifdef MINIMP3_FLOAT_OUTPUT
        typedef float W_D_TYPE;
        #define DATA_SIZE sizeof(float)
        #define write_wave write_float_wav
    #else
        typedef short W_D_TYPE;
        #define DATA_SIZE  sizeof(short)
        #define write_wave write_float_wav
    #endif

   
    typedef struct {
        size_t num_samples;
        size_t channels;
        W_D_TYPE *samples;
        float sample_rate;
        long file_size;
    } audio_data;
 
    float print_ad(audio_data *data);
    void free_audio(audio_data *audio);
    void read_file(const char *filename, uint64_t *size, uint8_t **data);
    audio_data read_wav(const char *filename,long file_size);
    audio_data read_mp3(const char *filename,long file_size);
    audio_data auto_detect(const char *filename);

#endif 
 