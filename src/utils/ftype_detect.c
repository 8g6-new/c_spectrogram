#include "../../headers/utils/ftype_detect.h"

static const char* mime_types_map[AUDIO_TYPE_COUNT] = {
    "audio/unknown",
    "audio/mpeg",
    "audio/wav",
    "audio/flac",
    "audio/ogg",
    "audio/aac",
    "audio/amr",
    "audio/opus"
};

audio_type detect_audio_type(const char *filename) {
    


    FILE *file = fopen(filename, "rb");
    if (!file) {
        perror("Error opening file");
        return AUDIO_UNKNOWN;
    }

    uint8_t buffer[MAX_HEADER_SIZE];
    size_t read_bytes = fread(buffer, 1, MAX_HEADER_SIZE, file);
    fclose(file);

    if (read_bytes < 12) {
        return AUDIO_UNKNOWN;
    }

    uint64_t header64;
    uint32_t header32;
    memcpy(&header64, buffer, sizeof(header64));
    memcpy(&header32, buffer, sizeof(header32));

    int is_wav   = (header32 == 0x46464952) & (*(uint32_t*)(buffer + 8) == 0x45564157); // WAV: "RIFF" + "WAVE"
    int is_mp3   = ((header32 & 0xFFFFFF) == 0x334449) | (((buffer[0] & 0xFF) == 0xFF) & ((buffer[1] & 0xE0) == 0xE0)); // MP3: "ID3" or MPEG frame
    int is_flac  = (header32 == 0x43614C66); // FLAC: "fLaC"
    int is_ogg   = (header32 == 0x5367674F); // OGG: "OggS"
    int is_opus  = is_ogg & (buffer[28] == 1); // OPUS: Ogg with Opus codec flag
    int is_aac   = ((buffer[0] == 0xFF) & ((buffer[1] & 0xF0) == 0xF0)); // AAC: ADTS Sync word
    int is_amr   = (*(uint32_t*)buffer == 0x524D4123); // AMR: "#!AMR"


    audio_type type =  is_wav * AUDIO_WAV |
           is_mp3 * AUDIO_MPEG |
           is_flac * AUDIO_FLAC |
           is_opus * AUDIO_OPUS |
           is_ogg * AUDIO_OGG |
           is_aac * AUDIO_AAC |
           is_amr * AUDIO_AMR |
           AUDIO_UNKNOWN;

    printf("%s auto detected to be %s\n",filename,mime_types_map[type]);

    return type;
}

inline const char* get_mime_type(audio_type type) {
    
    return mime_types_map[type];
}

