#ifndef FTYPE_DETECT_H
#define FTYPE_DETECT_H

    #include <stdio.h>
    #include <stdint.h>
    #include <string.h>
  
   
    #define MAX_HEADER_SIZE 36

    typedef enum {
        AUDIO_UNKNOWN=0,
        AUDIO_WAV,
        AUDIO_MPEG, 
        AUDIO_FLAC,
        AUDIO_OPUS,
        AUDIO_OGG,
        AUDIO_AAC,
        AUDIO_AMR,
        AUDIO_TYPE_COUNT
    } audio_type;


    extern const char* mime_types_map[AUDIO_TYPE_COUNT];
    
    audio_type detect_audio_type(const char *filename,long *file_size);
    const char* get_mime_type(audio_type type);

#endif 