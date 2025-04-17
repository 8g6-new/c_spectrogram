#include "../../headers/utils/ftype_detect.h"

const char* mime_types_map[AUDIO_TYPE_COUNT] = {
    "audio/unknown",  // AUDIO_UNKNOWN
    "audio/wav",      // AUDIO_WAV
    "audio/mpeg",     // AUDIO_MPEG 
    "audio/flac",     // AUDIO_FLAC
    "audio/opus",     // AUDIO_OPUS 
    "audio/ogg",      // AUDIO_OGG
    "audio/aac",      // AUDIO_AAC
    "audio/amr"       // AUDIO_AMR
};

audio_type detect_audio_type(const char *filename, long *file_size) {
    if (!filename) {
        return AUDIO_UNKNOWN;
    }

    FILE *file = fopen(filename, "rb");
    if (!file) {
        perror("Error opening file");
        return AUDIO_UNKNOWN;
    }

    if (fseek(file, 0, SEEK_END) != 0) {
        perror("Error seeking file end");
        fclose(file);
        return AUDIO_UNKNOWN;
    }

    *file_size = ftell(file);
    if (*file_size < 0) {
        perror("Error getting file size");
        fclose(file);
        return AUDIO_UNKNOWN;
    }

    if (fseek(file, 0, SEEK_SET) != 0) {
        perror("Error seeking file beginning");
        fclose(file);
        return AUDIO_UNKNOWN;
    }

    uint8_t buffer[MAX_HEADER_SIZE];
    size_t read_bytes = fread(buffer, 1, MAX_HEADER_SIZE, file);
    fclose(file);

    if (read_bytes < 12) {
        return AUDIO_UNKNOWN;
    }

    uint32_t riff_sig, wave_sig, flac_sig, ogg_sig, amr_sig;
    
    memcpy(&riff_sig, buffer, 4);
    memcpy(&wave_sig, buffer + 8, 4);
    memcpy(&flac_sig, buffer, 4);
    memcpy(&ogg_sig, buffer, 4);
    memcpy(&amr_sig, buffer, 4);


    // Check for WAV: "RIFF" + "WAVE"
    int is_wav  = (riff_sig == 0x46464952) && (wave_sig == 0x45564157);
    // Check for MP3: "ID3" or MPEG frame
    int is_mp3  = ((buffer[0] == 0x49 && buffer[1] == 0x44 && buffer[2] == 0x33) || ((buffer[0] & 0xFF) == 0xFF && (buffer[1] & 0xE0) == 0xE0));
    // Check for FLAC: "fLaC"
    int is_flac = (flac_sig == 0x43614C66);
    // Check for OGG: "OggS"
    int is_ogg  = (ogg_sig == 0x5367674F);
    // Check for OPUS: "OggS" followed by "OpusHead"
    int is_opus = is_ogg && (memcmp(buffer + 28, "OpusHead", 8) == 0);
    // Check for AAC: ADTS Sync word
    int is_aac  = ((buffer[0] == 0xFF) && ((buffer[1] & 0xF0) == 0xF0));
    // Check for AMR: "#!AMR"
    int is_amr  = (amr_sig == 0x524D4123);

    audio_type type = AUDIO_UNKNOWN;

               type =   is_wav *  AUDIO_WAV   |
                        is_mp3 *  AUDIO_MPEG  |
                        is_flac * AUDIO_FLAC  |
                        is_opus * AUDIO_OPUS  |
                        is_ogg *  AUDIO_OGG   |
                        is_aac *  AUDIO_AAC   |
                        is_amr *  AUDIO_AMR;

    printf("%s auto detected to be %s\n", filename, mime_types_map[type]);

    return type;
}


inline const char* get_mime_type(audio_type type) {
    if (type >= 0 && type < AUDIO_TYPE_COUNT) {
        return mime_types_map[type];
    }
    return mime_types_map[0]; 
}