#ifndef FTYPE_DETECT_H
#define FTYPE_DETECT_H

    #include <stdio.h>
    #include <stdint.h>
    #include <string.h>
  
   
    #define MAX_HEADER_SIZE 16

    typedef enum {
        AUDIO_UNKNOWN = 0,
        AUDIO_MPEG,
        AUDIO_WAV,
        AUDIO_FLAC,
        AUDIO_OGG,
        AUDIO_AAC,
        AUDIO_AMR,
        AUDIO_OPUS,
        AUDIO_TYPE_COUNT
    } audio_type;

    static const char* mime_types_map[AUDIO_TYPE_COUNT];
    audio_type detect_audio_type(const char *filename);
    const char* get_mime_type(audio_type type);

#endif 