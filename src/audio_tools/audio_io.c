#include "../../headers/audio_tools/audio_io.h"
#include "../../headers/audio_tools/minimp3.h"
#include "../../headers/utils/bench.h"

/*
 * The MIT License (MIT)
 * 
 * Copyright © 2025 Devadut S Balan
 * 
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the “Software”), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 * 
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */


float print_ad(audio_data *data) {
    if (!data || !data->samples) {
        printf("{\"error\": \"Invalid audio data\"}\n");
        return 0.0f;
    }
    
    float duration = (float)data->num_samples / ((float)(data->sample_rate * data->channels));
    
    printf("duration:%.3f\n",       duration);
    printf("channels:%zu\n",        data->channels);
    printf("num_samples:%zu\n",     data->num_samples);
    printf("sample_rate:%.4f\n",    data->sample_rate);
    printf("file_size_bytes:%ld\n", data->file_size);

    return duration;
}

void free_audio(audio_data *audio) {
    if (audio) {
        if (audio->samples) {  // Check if samples is not NULL before freeing
            free(audio->samples);
            audio->samples = NULL;
        }
    }
}

void read_file(const char *filename, uint64_t *size, uint8_t **data) {
    *size = 0;
    *data = NULL;

    FILE *fin = fopen(filename, "rb");
    
    if (!fin) {
        printf("\nError opening input file\n");
        perror("fopen");
        return;
    }

    struct stat st;
    if (fstat(fileno(fin), &st) != 0) {
        perror("fstat");
        fclose(fin);
        return;
    }

    *size = (uint64_t)st.st_size;

    *data = malloc(*size);

    if (!*data) {
        fprintf(stderr, "Memory allocation failed\n");
        *size = 0;
        fclose(fin);
        return;
    }

    if (fread(*data, 1, *size, fin) != *size) {
        perror("fread");
        free(*data);
        *data = NULL;
        *size = 0;
    }

    fclose(fin);
}

audio_data read_wav(const char *filename, long file_size) {
    audio_data audio = {0};
    audio.file_size = file_size;

    SNDFILE *file;
    SF_INFO sf_info = {0};  // Initialize sf_info to zero

    file = sf_open(filename, SFM_READ, &sf_info);
    if (!file) {
        fprintf(stderr, "Error opening file\n");
        return audio;
    }

    audio.num_samples = (size_t)sf_info.frames * sf_info.channels;
    audio.samples = (float*)malloc(audio.num_samples * sizeof(float));

    if (!audio.samples) {
        fprintf(stderr, "Memory allocation failed\n");
        sf_close(file);
        return audio;
    }

    sf_count_t frames_read = sf_readf_float(file, audio.samples, sf_info.frames);
    if (frames_read < sf_info.frames) {
        fprintf(stderr, "Error reading audio data: read %lld of %lld frames\n", 
                (long long)frames_read, (long long)sf_info.frames);
        // Adjust num_samples based on what was actually read
        audio.num_samples = (size_t)frames_read * sf_info.channels;
    }

    audio.sample_rate = sf_info.samplerate;
    audio.channels = sf_info.channels;

    sf_close(file);
    return audio;
}

audio_data read_mp3(const char *filename, long file_size) {
    audio_data audio = {0};
    audio.file_size = file_size;

    static mp3dec_t mp3d;
    mp3dec_init(&mp3d);

    uint64_t buf_size = 0;
    uint8_t *input_buf = NULL;
    
    START_TIMING();
    read_file(filename, &buf_size, &input_buf);
    END_TIMING("file_read");
   
    START_TIMING();
    if (!input_buf) {
        fprintf(stderr, "Failed to read input file: %s\n", filename);
        return audio;
    }

    // More conservative estimation of max samples to prevent buffer overflow
    size_t max_pcm_samples = (buf_size * MINIMP3_MAX_SAMPLES_PER_FRAME) / 128;
    size_t pcm_bsiz = max_pcm_samples * sizeof(W_D_TYPE) * 2;  // Account for stereo

    audio.samples = malloc(pcm_bsiz);
    if (!audio.samples) {
        fprintf(stderr, "Memory allocation failed\n");
        free(input_buf);
        return audio;
    }

    W_D_TYPE *full_pcm = (W_D_TYPE *)audio.samples;

    uint64_t decoded_samples = 0;
    int remaining_size = buf_size;
    uint8_t *input_ptr = input_buf;

    while (remaining_size > 0) {
        mp3dec_frame_info_t info;
        W_D_TYPE pcm[MINIMP3_MAX_SAMPLES_PER_FRAME * 2];  // Buffer for stereo

        int samples = mp3dec_decode_frame(&mp3d, input_ptr, remaining_size, pcm, &info);

        if (info.frame_bytes == 0 || remaining_size < info.frame_bytes) {
            break;
        }

        if (samples > 0) {
            // Check for potential buffer overflow
            if ((decoded_samples + samples) * info.channels > max_pcm_samples) {
                fprintf(stderr, "PCM buffer overflow prevented\n");
                break;
            }

            size_t copy_size = samples * sizeof(W_D_TYPE) * info.channels;
            memcpy(full_pcm + (decoded_samples * info.channels), pcm, copy_size);
            decoded_samples += samples;

            audio.channels = info.channels;
            audio.sample_rate = info.hz;
        }

        input_ptr += info.frame_bytes;
        remaining_size -= info.frame_bytes;
    }

    // Set the correct number of samples
    audio.num_samples = decoded_samples * audio.channels;
    
    free(input_buf);
    END_TIMING("dec_mp3");
    return audio;
}

audio_data auto_detect(const char *filename) {
    audio_data audio = {0};
    
    START_TIMING();
    audio_type type = detect_audio_type(filename, &audio.file_size);
    END_TIMING("auto_det");

    switch (type) {
        case AUDIO_WAV:
            audio = read_wav(filename, audio.file_size);
            END_TIMING("dec_wav");
            break;
        case AUDIO_MPEG:
            audio = read_mp3(filename, audio.file_size);
            break;
        default:
            fprintf(stderr, "Unknown or unsupported audio format\n");
            break;
    }
    
    return audio;
}
