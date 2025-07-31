#ifndef AUDIO_IO_H
    #define AUDIO_IO_H

    #include <stdio.h>       
    #include <stdlib.h>    
    #include <stdint.h>      
    #include <sys/stat.h>    
    #include <sndfile.h>
    #include <time.h>
    #include <math.h>

    /**
     * @file audio_io.h
     * @brief Audio file I/O and decoding support (WAV, MP3)
     */

    /**
     * @brief MP3 Frame Size Calculation Notes
     *
     * The size of an MP3 frame in bytes is given by the formula:
     *   frame_size = (K * bitrate) / sample_rate + padding
     *
     * K values:
     *   - MPEG-1 Layer III    → 144
     *   - MPEG-2/2.5 Layer III→ 72
     *
     * Worst-case estimation:
     *   - bitrate     = 32000 bps
     *   - sample_rate = 22050 Hz
     *   - padding     = 1
     *   - frame_size ≈ 209 bytes
     */

    #define WORST_BITRATE_BPS     32000
    #define WORST_SAMPLE_RATE     22050
    #define WORST_PADDING         1
    #define LAYER_III_FRAME_CONST 144

    /**
     * @brief Worst-case MP3 frame size (~209 bytes)
     */
    #define WORST_CASE_FRAME_SIZE ((LAYER_III_FRAME_CONST * WORST_BITRATE_BPS) / WORST_SAMPLE_RATE + WORST_PADDING)

    #include "../utils/ftype_detect.h"

    #define MINIMP3_ONLY_MP3
    #define MINIMP3_USE_SIMD
    #define MINIMP3_IMPLEMENTATION
    #define write_wave write_float_wav

    #ifdef MINIMP3_FLOAT_OUTPUT
        typedef float PARAM_DATATYPE;
    #else
        typedef short PARAM_DATATYPE;
    #endif

    /**
     * @brief Raw file buffer
     */
    typedef struct {
        uint8_t *data;   /**< Pointer to raw bytes */
        uint64_t size;   /**< Buffer size in bytes */
    } file_buffer;

    /**
     * @brief MP3 frame metadata
     */
    typedef struct {
        unsigned short *data;       /**< Frame deltas (packed) */
        int count;                  /**< Number of MP3 frames */
        float avg_byte_per_frame;  /**< Average frame size */
    } frames;

    /**
     * @brief Decoded audio structure
     */
    typedef struct {
        size_t num_samples;         /**< Total number of samples */
        size_t channels;            /**< Number of audio channels */
        PARAM_DATATYPE *samples;    /**< Interleaved sample data */
        float sample_rate;          /**< Sampling rate */
        long file_size;             /**< File size in bytes */
    } audio_data;

    /**
     * @brief Print audio data summary (e.g. stats)
     * @param data Pointer to audio_data
     * @return float Debug metric or duration (user-defined)
     */
    float print_ad(audio_data *data);

    /**
     * @brief Free allocated memory in audio_data
     * @param audio Pointer to audio_data
     */
    void free_audio(audio_data *audio);

    /**
     * @brief Read file into memory buffer
     * @param filename Path to file
     * @return file_buffer Structure with file contents
     */
    file_buffer read_file(const char *filename);

    /**
     * @brief Read WAV file using libsndfile
     * @param filename Path to WAV file
     * @param file_size File size in bytes
     * @return audio_data Struct with decoded PCM samples
     */
    audio_data read_wav(const char *filename, long file_size);

    /**
     * @brief Read MP3 file using minimp3
     * @param filename Path to MP3 file
     * @param file_size File size in bytes
     * @return audio_data Struct with decoded PCM samples
     */
    audio_data read_mp3(const char *filename, long file_size);

    /**
     * @brief Auto-detect file type (WAV/MP3) and decode
     * @param filename Path to audio file
     * @return audio_data Struct with decoded audio
     */
    audio_data auto_detect(const char *filename);

    /**
     * @brief Locate MP3 frame start positions for future multi threding
     * @param buf Pointer to file_buffer containing MP3 data
     * @return frames Struct containing frame offsets and metadata
     */
    frames find_mp3_frame_offsets(file_buffer *buf);

#endif // AUDIO_IO_H
