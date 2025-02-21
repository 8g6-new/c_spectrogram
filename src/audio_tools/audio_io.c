#include "../../headers/audio_tools/audio_io.h"
#include "../../headers/audio_tools/minimp3.h"
    

void print_ad(audio_data *data) {
    printf("Audio Data:\n");
    printf(" Seconds : %f\n",data->num_samples / (data->channels * data->sample_rate));
    printf("  Number of Channels: %zu\n", data->channels);
    printf("  Sample Rate: %f Hz\n", data->sample_rate);
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


audio_data read_wav(const char *filename){

    audio_data audio = {0};
    
    SNDFILE *file;
    SF_INFO sf_info;

    file = sf_open(filename, SFM_READ, &sf_info);
    if (!file) {
        fprintf(stderr, "Error opening file\n");
        return audio ;
    }

    audio.num_samples = (size_t)sf_info.frames * sf_info.channels;
    audio.samples     = (float*)malloc(audio.num_samples * sizeof(float));

    if (!audio.samples) {
        fprintf(stderr, "Memory allocation failed\n");
        free(audio.samples);
        audio.samples = NULL;
        sf_close(file);
        return audio ;
    }

    if (sf_readf_float(file,audio.samples, sf_info.frames) < sf_info.frames) {
        fprintf(stderr, "Error reading audio data\n");
        free(audio.samples);
        audio.samples = NULL;
        sf_close(file);
        return audio ;
    }

    audio.sample_rate = sf_info.samplerate;
    audio.channels    = sf_info.channels;

    sf_close(file);

    return audio;
}

audio_data read_mp3(char *filename) {
    
    audio_data audio  = {0};

    static mp3dec_t mp3d;
    mp3dec_init(&mp3d);

    uint64_t buf_size  = 0;
    uint8_t *input_buf = NULL;

    read_file(filename, &buf_size, &input_buf);

    if (!input_buf) {
        fprintf(stderr, "Failed to read input file: %s\n", filename);
        return audio;
    }

    size_t data_size       = sizeof(W_D_TYPE);
    size_t max_pcm_samples = (buf_size * MINIMP3_MAX_SAMPLES_PER_FRAME) / 128; // (mp3 max size estimate)
    size_t pcm_bsiz        = max_pcm_samples * data_size * 2;

    audio.samples = malloc(pcm_bsiz);
    if (!audio.samples) {
        fprintf(stderr, "Memory allocation failed\n");
        free(input_buf);
        return audio;
    }

    W_D_TYPE *full_pcm = (W_D_TYPE *)audio.samples; 

    uint64_t decoded_samples = 0;
    size_t remaining_size    = buf_size;

    while (remaining_size > 0) {
        mp3dec_frame_info_t info;
        W_D_TYPE pcm[MINIMP3_MAX_SAMPLES_PER_FRAME * 2];

        int samples = mp3dec_decode_frame(&mp3d, input_buf, remaining_size, pcm, &info);

        if (info.frame_bytes == 0 || remaining_size < info.frame_bytes) {
            break;
        }

        if (samples > 0) {
            if ((decoded_samples + samples) * 2 > max_pcm_samples) {
                fprintf(stderr, "PCM buffer overflow prevented\n");
                break;
            }

            size_t copy_size = samples * data_size * info.channels;
            memcpy(full_pcm + (decoded_samples * info.channels), pcm, copy_size);

            audio.channels = info.channels;
            audio.sample_rate = info.hz;
            decoded_samples += samples;
        }

        input_buf += info.frame_bytes;
        remaining_size -= info.frame_bytes;
    }

    audio.num_samples = decoded_samples;

    return audio;
}


audio_data auto_detect(const char *filename){

    clock_t start = clock();

    audio_type type    =  detect_audio_type(filename);

    audio_data audio   = {0};

    switch (type) {
        case 1:
            audio = read_mp3(filename);
            break;
        case 2:
            audio = read_wav(filename);
            break;
        default:
            break;
    }

    clock_t end = clock();

    printf("%f ms for file detection and decoding\n",(float)(end - start) * pow(10,EXP) / CLOCKS_PER_SEC);

    return audio;
}