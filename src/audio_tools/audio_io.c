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
        if (audio->samples) { 
            free(audio->samples);
            audio->samples = NULL;
        }
    }
}

file_buffer read_file(const char *filename) {
    file_buffer result = {.data = NULL, .size = 0};

    FILE *fin = fopen(filename, "rb");
    if (!fin) {
        fprintf(stderr, "Error: Failed to open file '%s'\n", filename);
        return result;
    }

    struct stat st;
    if (fstat(fileno(fin), &st) != 0) {
        fprintf(stderr, "Error: fstat failed for file '%s'\n", filename);
        fclose(fin);
        return result;
    }

    result.size = (uint64_t)st.st_size;
    result.data = (uint8_t *)malloc(result.size);
    if (!result.data) {
        fprintf(stderr, "Error: Memory allocation failed for file '%s' (%lu bytes)\n", filename, result.size);
        fclose(fin);
        return result;
    }

    if (fread(result.data, 1, result.size, fin) != result.size) {
        fprintf(stderr, "Error: fread failed for file '%s'\n", filename);
        free(result.data);
        result.data = NULL;
        result.size = 0;
    }

    fclose(fin);
    return result;
}

audio_data read_wav(const char *filename, long file_size) {
    audio_data audio = {0};
    audio.file_size = file_size;

    SNDFILE *file;
    SF_INFO sf_info = {0};  

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
        audio.num_samples = (size_t)frames_read * sf_info.channels;
    }

    audio.sample_rate = sf_info.samplerate;
    audio.channels = sf_info.channels;

    sf_close(file);
    return audio;
}


frames find_mp3_frame_offsets(file_buffer *buf) {
    frames result = {.data = NULL, .count = 0, .avg_byte_per_frame = 0.0f};

    int offset = 0;
    int free_format_bytes = 0;
    int frame_bytes = 0;

  
    int max_frames = buf->size / WORST_CASE_FRAME_SIZE;   // seee audio_io.h line 5

    result.data = (unsigned short *)malloc(max_frames * sizeof(unsigned short));
    if (!result.data) {
        fprintf(stderr, "Memory allocation failed for %d frames\n", max_frames);
        return result;
    }

    int frame_index = 0;
    int total_bytes = 0;

    while (offset < buf->size && frame_index < max_frames) {
        int next_offset = mp3d_find_frame(buf->data + offset, buf->size - offset, &free_format_bytes, &frame_bytes);

        if (frame_bytes == 0 || next_offset < 0)
            break;

        int total_advance = next_offset + frame_bytes;

        result.data[frame_index++] = (unsigned short)total_advance;
        total_bytes += total_advance;
        offset += total_advance;
    }

    result.count = frame_index;
    if (frame_index > 0)
        result.avg_byte_per_frame = (float)total_bytes / frame_index;

    return result;
}



audio_data read_mp3(const char *filename, long file_size) {
    audio_data audio = {0};
    audio.file_size = file_size;

    static mp3dec_t mp3d;
    mp3dec_init(&mp3d);

    START_TIMING();
    file_buffer buf = read_file(filename);
    END_TIMING("file_read");

    if (!buf.data || buf.size == 0) {
        fprintf(stderr, "Failed to read input file: %s\n", filename);
        return audio;
    }

    START_TIMING();
    frames f = find_mp3_frame_offsets(&buf);
    END_TIMING("frame_scan");

    printf("Total frames: %d\n", f.count);
    printf("Average frame size: %.2f bytes\n", f.avg_byte_per_frame);


    // Estimate max PCM samples (safe upper bound at 32kbps)
    size_t max_pcm_samples = (buf.size * MINIMP3_MAX_SAMPLES_PER_FRAME) / 32;
    size_t pcm_bsiz        = max_pcm_samples * sizeof(PARAM_DATATYPE) * 2;

    audio.samples = malloc(pcm_bsiz);
    if (!audio.samples) {
        fprintf(stderr, "Memory allocation failed\n");
        free(buf.data);
        free(f.data);
        return audio;
    }

    PARAM_DATATYPE *full_pcm = (PARAM_DATATYPE *)audio.samples;
    uint64_t decoded_samples = 0;

    uint8_t *input_ptr = buf.data;
    int remaining_size = buf.size;

    START_TIMING();
    while (remaining_size > 0) {
        mp3dec_frame_info_t info;
        PARAM_DATATYPE pcm[MINIMP3_MAX_SAMPLES_PER_FRAME * 2];

        int samples = mp3dec_decode_frame(&mp3d, input_ptr, remaining_size, pcm, &info);

        if (info.frame_bytes == 0 || remaining_size < info.frame_bytes)
            break;

        if (samples > 0) {
            if ((decoded_samples + samples) * info.channels > max_pcm_samples) {
                fprintf(stderr, "PCM buffer overflow prevented\n");
                break;
            }

            size_t copy_size = samples * sizeof(PARAM_DATATYPE) * info.channels;
            memcpy(full_pcm + (decoded_samples * info.channels), pcm, copy_size);
            decoded_samples += samples;

            audio.channels = info.channels;
            audio.sample_rate = info.hz;
        }

        input_ptr += info.frame_bytes;
        remaining_size -= info.frame_bytes;
    }
    END_TIMING("dec_mp3");

    audio.num_samples = decoded_samples * audio.channels;

    free(buf.data);
    free(f.data);

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
