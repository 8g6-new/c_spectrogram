#ifndef AUDIO_TOOLS_H
#define AUDIO_TOOLS_H
    #include <sndfile.h>
    #include <fftw3.h>
    #include <stddef.h>
    #include <stdbool.h>


    #define EXP 6 // Microseconds
    #define MAX_FILTERS 128 // Max filters for mel spectrogram

    typedef struct {
        size_t num_samples;
        size_t channels;
        size_t frames;
        float *samples;
        float sample_rate;
    } audio_data;

    typedef struct {
        size_t output_size;
        size_t num_frequencies;
        float *magnitudes;
        float *phases;
        float fft_times;
        float *phasers;
    } stft_d;
    

    audio_data read_file(const char *filename);

    static void init_fft_output(stft_d *result, unsigned int window_size, unsigned int hop_size, unsigned int num_samples);
    stft_d stft(float *samples, size_t window_size, size_t hop_size, uint64_t num_samples, float *window_values);
    void free_stft(stft_d *result);


    void spectrogram(stft_d *result, const char *output_file, unsigned char bg_clr[4], bool db);
    void mel_spectrogram(stft_d *result,const char *output_file,unsigned char bg_clr[4],bool db,size_t num_filters,float *mel_filter_bank,float *mel_vales);
    void MFCC(float *mel_values,const char *output_file,size_t w,unsigned char bg_clr[4],bool db,size_t num_filters,size_t num_coff);

    void window_function(float *window_values, size_t window_size, const char *window_type);
    void mel_filter(float min_f, float max_f, size_t n_filters, float sr, size_t fft_size, float *filter);

    static double hz_to_mel(double f, float mid_f);
    static double mel_to_hz(double m, float mid_f);
    void calculate_frequencies(float *frequencies, size_t window_size, float sample_rate);
    size_t safe_diff(size_t a, size_t b);

#endif