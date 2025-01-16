#include "audio_io.h"

void readFile(const char *filename, uint64_t *size, uint8_t **data){

}



audio_data wav_to_samples(const char *filename){

}

audio_data mp3_to_samples(char *input_filename){

}

int ends_with_wav_or_mp3(const char *filename) {

    const char *end = filename;
    while (*end) ++end; 

    if (end - filename >= 4 && 
        end[-1] == '3' && end[-2] == 'p' && end[-3] == 'm' && end[-4] == '.') {
        return 1; 
    }

    if (end - filename >= 4 && 
        end[-1] == 'v' && end[-2] == 'a' && end[-3] == 'w' && end[-4] == '.') {
        return 1; 
    }

    return 0;
}
    
audio_data read_audio_file(const char *filename){

}

