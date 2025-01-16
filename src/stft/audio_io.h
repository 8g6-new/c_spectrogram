#ifndef AUDIO_IO_H
#define AUDIO_IO_H

#include <stdint.h>
#include "minimp3.h"

// Supported audio file extensions
    char exts[2][5] = {
        ".mp3",
        ".wav"
    };

    // Define audio data type based on output format
    #define MINIMP3_ONLY_MP3
    #define MINIMP3_USE_SIMD
    #define MINIMP3_IMPLEMENTATION
    #define MINIMP3_FLOAT_OUTPUT

    #ifdef MINIMP3_FLOAT_OUTPUT
        typedef float wave_data_type;
    #else
        typedef short wave_data_type;
    #endif

    typedef struct {
        unsigned int sample_rate;     /* Sample rate of the audio */
        unsigned int channels;        /* Number of channels */
        long frames;                  /* Total number of frames in the file (samples per channel) */
        wave_data_type *samples;      /* Audio data */
    } audio_data;

    audio_data wav_to_samples(const char *filename);
    audio_data mp3_to_samples(char *input_filename);
    void readFile(const char *filename, uint64_t *size, uint8_t **data);
    const short ends_with(const char *filename);
    audio_data read_audio_file(const char *filename);

#endif // AUDIO_IO_H
