/**
 * @file spectral_features.h
 * @brief Spectral feature extraction (STFT, windowing, filter banks).
 * 
 * This module provides core routines for computing short-time Fourier transforms,
 * generating and applying filter banks (e.g., Mel, Bark, ERB), and converting frequency
 * scales used in auditory signal processing.
 * 
 * @ingroup audio_features
 */

#ifndef SPECTRAL_FEATURES_H
    #define SPECTRAL_FEATURES_H

    #include <stddef.h>
    #include <stdbool.h>
    #include <math.h>
    #include <fftw3.h>
    #include <string.h>
    #include <omp.h>
    #include <cblas.h>

    #include "audio_io.h"
    #include "../utils/bench.h"


     #define LOG10_E 0.4342944819

    /** @addtogroup audio_features
    *  @{
    */

    /**
    * @brief Supported filter bank types
    */
    typedef enum {
        F_MEL,       /**< Mel scale */
        F_BARK,      /**< Bark scale */
        F_ERB,       /**< Equivalent Rectangular Bandwidth */
        F_CHIRP,     /**< Chirp scale */
        F_CAM,       /**< Cambridge ERB-rate */
        F_LOG10,     /**< Log10 frequency scaling */
        F_UNKNOWN    /**< Invalid/unsupported filter type */
    } filter_type_t;

    static const char *FILTER_TYPE_NAMES[] = {
        [F_MEL]    = "MEL",
        [F_BARK]   = "BARK",
        [F_ERB]    = "ERB",
        [F_CHIRP]  = "CHIRP",
        [F_CAM]    = "CAM",
        [F_LOG10]  = "LOG10",
        [F_UNKNOWN]= "UNKNOWN"
    };


    /**
    * @brief STFT data structure
    */
    typedef struct {
        float sample_rate;
        size_t output_size;
        size_t num_frequencies;
        size_t total_length;
        float *frequencies;
        float *magnitudes;
        float *phasers;
    } stft_d;

    /**
    * @brief FFT plan wrapper
    */
    typedef struct {
        fftwf_plan plan;
        fftwf_complex *in;
        fftwf_complex *out;
    } fft_d;

    /**
    * @brief Filter bank representation
    */
    typedef struct {
        size_t  *freq_indexs;
        float   *weights;
        size_t   size;
        const size_t num_filters;
    } filter_bank_t;

    /* ---- STFT and FFT Utilities ---- */

    bool init_fft_output(stft_d *result, unsigned int window_size, unsigned int hop_size, unsigned int num_samples);
    fft_d init_fftw_plan(const size_t window_size, const char *cache_dir);
    stft_d stft(audio_data *audio, const size_t window_size, const size_t hop_size, float *window_values, fft_d *fft);
    void free_stft(stft_d *result);
    void free_fft_plan(fft_d *fft);
    void cleanup_fft_threads(fft_d *thread_ffts, const size_t num_threads);
    void calculate_frequencies(stft_d *result, size_t window_size, float sample_rate);

    /* ---- Window and Filterbank Utilities ---- */

    void window_function(float *window_values, size_t window_size, const char *window_type);
    filter_bank_t gen_filterbank(filter_type_t type, float min_f, float max_f, size_t n_filters, float sr, size_t fft_size, float *filter);
    filter_type_t parse_filter_type(const char *name);
    void print_melbank(const filter_bank_t *v);

    /* ---- Frequency and Scale Utilities ---- */

    size_t hz_to_index(size_t num_freq, size_t sample_rate, float f);
    float safe_diff(size_t a, size_t b);
    float brachless_db(float mag, bool db);

    /* ---- Frequency Scale Conversion ---- */

    double hz_to_mel(double hz);     double mel_to_hz(double mel);
    double hz_to_bark(double hz);    double bark_to_hz(double bark);
    double hz_to_erb(double hz);     double erb_to_hz(double erb);
    double hz_to_chirp(double hz);   double chirp_to_hz(double chirp);
    double hz_to_cam(double hz);     double cam_to_hz(double cam);
    double hz_to_log10(double hz);   double log10_to_hz(double val);
    double hz_to_cent(double hz);    double cent_to_hz(double cent);

    /** @} */  // end of audio_features

#endif // SPECTRAL_FEATURES_H
