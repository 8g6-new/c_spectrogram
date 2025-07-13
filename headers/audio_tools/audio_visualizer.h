#ifndef AUDIO_VISUALIZER_H
#define AUDIO_VISUALIZER_H 
    
    #include "audio_io.h"
    #include "../libheatmap/heatmap_tools.h"
    #include "spectral_features.h"
    
    typedef struct {
        float start_f;
        float end_f;
        size_t start_d;
        size_t end_d;
    } bounds_t;
    
    typedef struct{
        bounds_t time;
        bounds_t freq;
    } bounds2d_t;

    typedef struct{
        char bg_color[4];
        const unsigned short int cs_enum;
        const bool db;
        const char output_file[500];
    } plot_t;


    typedef struct{
        float *coeffs;
        const size_t num_filters;
        const size_t num_coff;
    } mffc_t;


    void   print_bounds(bounds2d_t *bounds);
    void   init_bounds(bounds2d_t *bounds,stft_d *result);
    void   fast_copy(float *contious_mem,float *mags,bounds2d_t *bounds,const size_t length);
    void   set_limits(bounds2d_t *bounds,const size_t max_freq,const size_t max_time);
    void   spectrogram(float *contious_mem,bounds2d_t *bounds,plot_t *settings);

    float *mel_spectrogram(float *contious_mem,const size_t num_filters,const size_t num_freq,float *mel_filter_bank, bounds2d_t *bounds, plot_t *settings);
    mffc_t precompute_cosine_coeffs(const size_t num_filters,const size_t num_coff);
    void   mfcc(float *mel_values, mffc_t *dft_coff,bounds2d_t *bounds, plot_t *settings);

#endif


