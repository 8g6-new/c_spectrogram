#ifndef AUDIO_VISUALIZER_H
    #define AUDIO_VISUALIZER_H

    #include <omp.h>

    #include "audio_io.h"
    #include "../libheatmap/heatmap_tools.h"
    #include "spectral_features.h"

    /**
     * @brief Frequency/time bound range
     */
    typedef struct {
        float start_f;     /**< Start frequency */
        float end_f;       /**< End frequency */
        size_t start_d;    /**< Start index (time or frequency) */
        size_t end_d;      /**< End index (time or frequency) */
    } bounds_t;

    /**
     * @brief 2D time-frequency bounds
     */
    typedef struct {
        bounds_t time;     /**< Time bounds */
        bounds_t freq;     /**< Frequency bounds */
    } bounds2d_t;

    /**
     * @brief Plotting configuration
     */
    typedef struct {
        char bg_color[4];              /**< Background RGBA */
        const unsigned short int cs_enum; /**< Colormap selection enum */
        const bool db;                /**< Whether to use dB scaling */
        char output_file[500]; /**< Output file name */
        size_t w;                     /**< Width */
        size_t h;                     /**< Height */
    } plot_t;

    /**
     * @brief DCT coefficient structure
     */
    typedef struct {
        float *coeffs;            /**< Cosine coefficients */
        const size_t num_filters; /**< Number of Mel filters */
        const size_t num_coff;    /**< Number of DCT coeffs */
    } dct_t;

    /**
     * @brief Print the bounds for debugging
     * @param bounds Pointer to bounds2d_t
     */
    void print_bounds(bounds2d_t *bounds);

    /**
     * @brief Set maximum limits for bounds
     * @param bounds Pointer to bounds2d_t
     * @param max_freq Maximum frequency index
     * @param max_time Maximum time index
     */
    void set_limits(bounds2d_t *bounds, const size_t max_freq, const size_t max_time);

    /**
     * @brief Initialize bounds using STFT result
     * @param bounds Pointer to bounds2d_t
     * @param result Pointer to stft_d
     */
    void init_bounds(bounds2d_t *bounds, stft_d *result);

    /**
     * @brief Copy a subregion from data to mags based on bounds
     * @param data Input data array
     * @param mags Output magnitude array
     * @param bounds Pointer to bounds2d_t
     * @param length Total length of input data
     */
    void fast_copy(float *data, float *mags, bounds2d_t *bounds, const size_t length);

    /**
     * @brief Plot a heatmap
     * @param data Data array to plot
     * @param bounds Pointer to bounds2d_t
     * @param settings Pointer to plot_t settings
     */
    void plot(float *data, bounds2d_t *bounds, plot_t *settings);

    /**
     * @brief Apply Mel/Bark filter bank to data
     * @param data Input FFT data
     * @param num_filters Number of filters
     * @param num_freq Number of frequency bins
     * @param mel_filter_bank Precomputed filter bank
     * @param bounds Pointer to bounds2d_t
     * @param settings Optional plot settings
     * @return float* Pointer to filter bank output (allocated)
     */
    float *apply_filter_bank(float *data, const size_t num_filters, const size_t num_freq, float *mel_filter_bank, bounds2d_t *bounds, plot_t *settings);

    /**
     * @brief Generate DCT cosine coefficients
     * @param num_filters Number of filters
     * @param num_coff Number of output DCT coefficients
     * @return dct_t Struct containing DCT coefficient array
     */
    dct_t gen_cosine_coeffs(const size_t num_filters, const size_t num_coff);

    /**
     * @brief Compute Filtered Cepstral Coefficients (FCC)
     * @param mel_values Input Mel spectrum
     * @param dft_coff Pointer to DCT coefficient struct
     * @param bounds Pointer to bounds2d_t
     * @param settings Optional plot settings
     * @return float* Array of cepstral coefficients
     */
    float *FCC(float *mel_values, dct_t *dft_coff, bounds2d_t *bounds, plot_t *settings);

#endif // AUDIO_VISUALIZER_H
