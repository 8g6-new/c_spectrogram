# c_spectrogram : Audio Signal Processing and Visualization


This repository provides a high-performance C library for audio signal processing and visualization, featuring Short-Time Fourier Transform (STFT), Mel spectrograms, Mel-Frequency Cepstral Coefficients (MFCC), and heatmap visualizations. Designed for large-scale audio datasets, it leverages optimizations like SIMD (SSE, SSE2, AVX, AVX2), FFTW with wisdom caching, OpenMP parallelization, and BLAS for fast matrix operations. The library supports WAV and MP3 inputs, produces professional-grade visualizations with customizable color schemes, and includes robust benchmarking tools for performance analysis.

## Features

- **Audio I/O**:
  - Supports WAV (via `libsndfile`) and MP3 (via `minimp3`) formats.
  - Automatic file type detection for seamless input handling.

- **Short-Time Fourier Transform (STFT)**:
  - Computes STFT with FFTW, optimized using wisdom caching for fast FFT planning.
  - Supports multiple window functions (Hann, Hamming, Blackman, Bartlett, Blackman-Harris, Flat-top, Gaussian, Kaiser).
  - SIMD optimizations (SSE, SSE2, AVX, AVX2) for magnitude and phase calculations.
  - Configurable window size, hop size, and frequency bounds.

- **Mel Spectrogram**:
  - Generates Mel spectrograms using dynamically computed Mel filter banks.
  - Optimized with BLAS (`cblas_sdot`) for fast filter bank application.
  - Parallelized with OpenMP for scalability on multi-core CPUs.
  - Optional dB scaling with branchless computation.

- **Mel-Frequency Cepstral Coefficients (MFCC)**:
  - Computes MFCCs via precomputed DCT coefficients and BLAS-accelerated operations.
  - Parallelized with OpenMP for efficient processing.
  - Visualizes MFCCs as heatmaps with customizable color schemes.

- **Visualization**:
  - Renders spectrograms, Mel spectrograms, and MFCCs as PNG heatmaps using a [libheatmap](https://github.com/lucasb-eyer/libheatmap).
  - Supports extensive color schemes (e.g., Blues, Viridis, Jet, Inferno) in discrete, mixed, and soft variants, with built-in and OpenCV-like options.
  - Configurable time/frequency bounds for focused visualizations.
  - Optimized memory copying for cache-friendly visualization.

- **Benchmarking**:
  - Detailed performance profiling with microsecond precision for key functions (e.g., STFT, Mel spectrogram, MFCC, plotting).
  - Ranked timing output with color-coded visualization and percentage of total runtime.
  - JSON and raw timing formats for integration with analysis tools.

- **Performance Optimizations**:
  - SIMD (SSE, SSE2, AVX, AVX2) with Fused Multiply-Add (FMA) and 16-bit floating-point conversions.
  - OpenMP parallelization for STFT, Mel spectrogram, MFCC, and visualization loops.
  - FFTW wisdom caching to reuse optimized FFT plans.
  - BLAS integration for fast linear algebra in Mel and MFCC computations.
  - Aggressive compiler optimizations (e.g., `-ffast-math`, `-march=native`, `-funroll-loops`, Link-Time Optimization).
  - Aligned memory allocations for SIMD and cache efficiency.

- **Applications**:
  - Ideal for bioacoustics (e.g., bird call analysis, tested with `black_woodpecker.wav`, `bird.mp3`).
  - Suitable for large-scale audio processing, feature extraction for machine learning, and DSP research.

## Performance Notes

- **Optimized Components**: MP3 decoding and plot saving outperform Python equivalents due to efficient I/O and PNG handling.
- **Scalability**: OpenMP parallelization scales with CPU cores, and SIMD optimizations leverage modern hardware (AVX2, FMA).
- **Benchmarking Insights**: Detailed timing reports (e.g., `print_bench_ranked`) help identify bottlenecks, guiding further optimization.

## Requirements

- **Compiler**: GCC or Clang with C11 support.
- **Libraries**:
  - FFTW3 (`libfftw3f`) for fast Fourier transforms.
  - libsndfile for WAV file handling.
  - OpenMP for parallel processing.
  - BLAS (e.g., OpenBLAS) for matrix operations.
  - libpng for PNG output.
- **Optional**: OpenCV for additional color scheme support (`opencv_like` build target).
- **Hardware**: Supports SSE, SSE2, AVX, and AVX2 instructions 

## Installation

1. **Install Dependencies** (Ubuntu/Debian example):
   ```bash
   sudo apt-get update
   sudo apt-get install libfftw3-dev libsndfile1-dev libopenblas-dev libpng-dev libomp-dev
   ```

2. **Clone the Repository**:
   ```bash
   git clone https://github.com/8g6-new/c_spectrogram
   cd c_spectrogram
   ```

3. **Build the Project**:
   - For built-in color schemes:
     ```bash
     make builtin
     ```
   - For OpenCV-like color schemes:
     ```bash
     make opencv_like
     ```
   - For a shared library:
     ```bash
     make shared
     ```
   - For debugging (built-in color schemes):
     ```bash
     make debug_builtin
     ```
   - For debugging (OpenCV-like color schemes):
     ```bash
     make debug_opencv_like
     ```

   The build process creates executables in the `build/builtin` or `build/opencv` directories and generates FFTW wisdom files in `cache/FFT` (e.g., `1024.wisdom`).

4. **Verify Build**:
   Check that the executable (`builtin` or `opencv_like`) is created and wisdom files are generated.

## Usage

### Command-Line Interface
The `main` program processes an audio file and generates STFT, Mel spectrogram, and MFCC visualizations with user-specified parameters.

```bash
./build/builtin/main <input_file> <output_prefix> <window_size> <hop_size> <window_type> <num_mel_banks> <min_mel> <max_mel> <num_mfcc_coeffs> <cs_stft> <cs_mel> <cs_mfcc> <cache_folder>
```

**Example**:
```bash
./build/builtin/main tests/files/black_woodpecker.wav outputs/black_woodpecker 2048 512 hann 40 20.0 8000.0 13 0 0 0 cache/FFT
```

This processes `black_woodpecker.wav` with:
- Window size: 2048 samples
- Hop size: 512 samples
- Hann window
- 40 Mel filters (20 Hz to 8 kHz)
- 13 MFCC coefficients
- Default color scheme (index 0)
- FFTW wisdom files in `cache/FFT`

**Output**:
- `outputs/black_woodpecker_stft.png` (STFT spectrogram)
- `outputs/black_woodpecker_mel.png` (Mel spectrogram)
- `outputs/black_woodpecker_mfcc.png` (MFCC heatmap)
- Console output with audio metadata (e.g., duration, sample rate) and ranked benchmark timings

### Code Example
Below is a minimal example of using the library programmatically:

```c
#include "headers/audio_tools/audio_io.h"
#include "headers/audio_tools/spectral_features.h"
#include "headers/audio_tools/audio_visualizer.h"
#include "headers/utils/bench.h"

int main() {
    // Initialize benchmarking
    benchmark_init();

    // Load audio
    audio_data audio = auto_detect("tests/files/black_woodpecker.wav");
    print_ad(&audio);

    // Compute STFT
    size_t window_size = 2048, hop_size = 512;
    float *window_values = malloc(window_size * sizeof(float));
    window_function(window_values, window_size, "hann");
    fft_d fft = init_fftw_plan(window_size, "cache/FFT");
    stft_d result = stft(&audio, window_size, hop_size, window_values, &fft);

    // Compute Mel filter bank
    size_t num_filters = 40;
    float *mel_filter_bank = calloc((result.num_frequencies + 1) * (num_filters + 2), sizeof(float));
    melbank_t melbank = mel_filter(20.0f, 8000.0f, num_filters, audio.sample_rate, window_size, mel_filter_bank);

    // Compute MFCC coefficients
    size_t num_coff = 13;
    mffc_t dft_coff = precompute_cosine_coeffs(num_filters, num_coff);

    // Prepare visualization bounds
    bounds2d_t bounds = { .freq = {20.0f, 8000.0f} };
    init_bounds(&bounds, &result);
    set_limits(&bounds, result.num_frequencies, result.output_size);

    // Copy STFT magnitudes
    float *contious_mem = malloc((bounds.freq.end_d - bounds.freq.start_d) * (bounds.time.end_d - bounds.time.start_d) * sizeof(float));
    fast_copy(contious_mem, result.magnitudes, &bounds, result.num_frequencies);

    // Visualize STFT spectrogram
    plot_t settings = { .db = true, .cs_enum = CS_Blues, .output_file = "outputs/stft.png" };
    settings.bg_color[0] = 0; settings.bg_color[1] = 0; settings.bg_color[2] = 0; settings.bg_color[3] = 255;
    spectrogram(contious_mem, &bounds, &settings);

    // Compute and visualize Mel spectrogram
    settings.output_file = "outputs/mel.png";
    float *mel_values = mel_spectrogram(contious_mem, num_filters, result.num_frequencies, mel_filter_bank, &bounds, &settings);

    // Compute and visualize MFCC
    settings.output_file = "outputs/mfcc.png";
    mfcc(mel_values, &dft_coff, &bounds, &settings);

    // Print benchmark results
    print_bench_ranked();

    // Cleanup
    free(mel_values);
    free(mel_filter_bank);
    free(contious_mem);
    free(window_values);
    free_stft(&result);
    free_fft(&fft);
    free_audio(&audio);
    free(dft_coff.coeffs);
    free(melbank.freq_indexs);
    free(melbank.weights);
    return 0;
}
```

## Sample Visualizations

Below are sample visualizations generated from `black_woodpecker.wav`, showcasing STFT spectrograms, Mel spectrograms, and MFCC heatmaps with various color schemes. These examples use the `builtin` and `opencv_like` builds.

### STFT Spectrograms (Built-in Color Schemes)
- **Blues (Soft)**: Smooth gradient for clear frequency visualization.
  ![STFT Blues Soft](outputs/colorschemes/libheatmap_defaults/soft/black_woodpecker_stft_Blues_soft.png)
- **Spectral (Discrete)**: High-contrast for distinct frequency bands.
  ![STFT Spectral Discrete](outputs/colorschemes/libheatmap_defaults/discrete/black_woodpecker_stft_Spectral_discrete.png)

### STFT Spectrograms (OpenCV-like Color Schemes)
- **Viridis**: Popular for scientific visualization, emphasizing frequency dynamics.
  ![STFT Viridis](outputs/colorschemes/opencv_like/images/black_woodpecker_stft_Viridis.png)
- **Jet**: Classic colormap for highlighting intensity variations.
  ![STFT Jet](outputs/colorschemes/opencv_like/images/black_woodpecker_stft_Jet.png)

### Mel Spectrogram and MFCC
- **Mel Spectrogram (Cividis)**: Visualizes Mel filter bank output.
  ![Mel Spectrogram](outputs/functions/black_woodpecke_mel.png)
- **MFCC (Blues Soft)**: Displays cepstral coefficients for feature extraction.
  ![MFCC](outputs/functions/black_woodpecke_mfcc.png)

To explore all available color schemes (e.g., Blues, Viridis, Jet, Inferno in discrete, mixed, mixed_exp, and soft variants), refer to the `README.MD` files in:
- [`outputs/colorschemes/libheatmap_defaults/README.MD`](./outputs/colorschemes/libheatmap_defaults/README.MD) for built-in color schemes.
- [`outputs/colorschemes/opencv_like/README.MD`](./outputs/colorschemes/opencv_like/README.MD) for OpenCV-like color schemes.

These files include comprehensive galleries of all color schemes applied to `black_woodpecker.wav`.

## 🎨 Colormap Enum Reference
All supported colormaps are listed in the file:

```bash
output/colors.json
```
This file maps human-readable names to internal enum IDs for both:

OpenCV-like colormaps (e.g., JET, VIRIDIS, HOT)

Built-in scientific colormaps (e.g., Blues.soft, Spectral.mixed_exp)

Refer [`outputs/README.MD`](./outputs/README.MD)


## Output Directory Structure
The `outputs` directory contains:
- `colorschemes/libheatmap_defaults`: STFT spectrograms with built-in color schemes in:
  - `discrete`: High-contrast, distinct color steps (e.g., `black_woodpecker_stft_Blues_discrete.png`).
  - `mixed`: Smooth transitions between colors (e.g., `black_woodpecker_stft_Blues_mixed.png`).
  - `mixed_exp`: Exponentially scaled mixed colors (e.g., `black_woodpecker_stft_Blues_mixed_exp.png`).
  - `soft`: Softened gradients for aesthetic visualization (e.g., `black_woodpecker_stft_Blues_soft.png`).
- `colorschemes/opencv_like/images`: STFT spectrograms with OpenCV-inspired color schemes (e.g., `black_woodpecker_stft_Viridis.png`, `black_woodpecker_stft_Jet.png`).
- `functions`: Mel spectrograms and MFCC heatmaps (e.g., `black_woodpecker_mel.png`, `black_woodpecker_mfcc.png`).

### Benchmarking Output
The `print_bench_ranked` function produces a ranked table of function execution times, including:
- Function name (e.g., `stft`, `mel`, `mfcc`, `stft_plot`).
- Execution time (in µs, ms, or s).
- Percentage of total runtime.
- Color-coded visualization for quick bottleneck identification.

Example output:
```
---------------------------------------------------------
| Function             | Exec Time    | % of total runtime |
---------------------------------------------------------
| stft                |   12.345 ms  |  45.6789% |
[▰▰▰▰▰▰▰▰▰▰▰▰▰▰▰▰▰▰▰▰        ]
| mel                 |    8.901 ms  |  32.1234% |
[▰▰▰▰▰▰▰▰▰▰▰▰▰▰                ]
| mfcc                |    3.456 ms  |  12.3456% |
[▰▰▰▰▰▰                      ]
...
---------------------------------------------------------
```

## Build Configuration

The `Makefile` provides flexible build targets with aggressive optimizations:

- **Compiler Flags**:
  - `-ffast-math`, `-march=native`, `-mtune=native` for CPU-specific optimizations.
  - `-funroll-loops`, `-floop-interchange`, `-floop-unroll-and-jam` for loop optimizations.
  - `-flto`, `-fuse-linker-plugin` for Link-Time Optimization.
  - `-mavx`, `-msse4.2`, `-mavx2`, `-mfma` for SIMD vectorization (AVX512 code commented out).
  - `-ftree-vectorize`, `-ftree-loop-vectorize` for automatic loop vectorization.

- **Debugging Flags**:
  - `-Og`, `-g`, `-fno-omit-frame-pointer` for debug builds.
  - `-fsanitize=address`, `-fsanitize=leak`, `-fsanitize=undefined` for memory and undefined behavior detection.

- **Build Targets**:
  - `make builtin`: Builds with built-in color schemes.
  - `make opencv_like`: Builds with OpenCV-like color schemes.
  - `make shared`: Builds a shared library (`libyourlib.so`).
  - `make debug_builtin` / `make debug_opencv_like`: Builds debug versions.
  - `make clean`: Removes build artifacts.

## Project Structure

```
.
├── cache/FFT/              # FFTW wisdom files for optimized FFT plans
├── headers/                # Header files for audio tools, heatmap, and utilities
├── outputs/                # Generated spectrograms, Mel spectrograms, and MFCC heatmaps
├── src/                    # Source code for audio processing, visualization, and utilities
│   ├── libheatmap/         # Heatmap visualization code
│   ├── png_tools/          # PNG output utilities
│   ├── utils/              # Benchmarking and utility functions
│   ├── audio_tools/        # Audio I/O, STFT, Mel, and MFCC computation
├── tests/files/            # Test audio files (e.g., black_woodpecker.wav, bird.mp3)
├── main.c                  # Main program for testing the pipeline
├── Makefile                # Build configuration
└── README.md               # This file
```

## Performance Optimization Tips

- **Multi-threading**: Adjust OpenMP threads with `OMP_NUM_THREADS` (e.g., `export OMP_NUM_THREADS=8`) to match your CPU cores.
- **SIMD**: Ensure `-mavx2` and `-mfma` are supported by your CPU. Uncomment AVX512 code in `spectral_features.c` for modern CPUs (e.g., Skylake-X).
- **BLAS**: Use an optimized BLAS implementation (e.g., OpenBLAS, MKL) for faster Mel and MFCC computations.
- **Wisdom Caching**: Pre-generate FFTW wisdom files for common window sizes to minimize planning overhead.
- **Matrix Operations**: Replace `cblas_sdot` loops in `mel_spectrogram` and `mfcc` with `cblas_sgemm` for batched matrix multiplications.
- **Fast Math Caution**: `-ffast-math` may introduce floating-point inaccuracies. Test thoroughly for numerical stability.
- **GPU Potential**: The pipeline is GPU-ready for STFT (cuFFT), Mel spectrograms (cuBLAS), and visualization (CUDA kernels). See "Future Work" for details.

## Future Work

- **GPU Acceleration**: Implement CUDA-based STFT (cuFFT), Mel spectrograms (cuBLAS), and visualization to match or exceed Librosa’s performance.
- **Sparse Matrix Operations**: Optimize Mel filter bank with sparse matrix formats (e.g., CSR) to reduce memory and computation.
- **Real-Time Processing**: Extend the pipeline for streaming audio analysis.
- **Enhanced Benchmarking**: Add memory usage and CPU/GPU utilization metrics to `bench.h`.
- **Simple Heatmap for Faster Plotting**: Replace the `libheatmap` `heatmap_add_weighted_point()` call used on every loop iteration, as it was identified as a significant bottleneck in image plotting. A custom, simplified heatmap generator will be implemented for faster rendering.
- **Memory Pooling and Allocation**: Implement memory pooling using a memory arena to improve memory utilization, prevent issues like use-after-free, and enhance overall memory management.
- **Advanced Memory Management**: Investigate and integrate state-of-the-art buddy allocators or other advanced memory management techniques to optimize memory allocation, reduce fragmentation, and improve performance in resource-constrained environments.
- **Documentation**: Add detailed API docs and usage examples to `headers/` and `README.md`.


## Acknowledgments

- Built with inspiration from [Librosa](https://librosa.org/), aiming for high-performance audio processing in C.

- Tested on bioacoustics datasets (e.g., bird calls), with thanks to open-source audio libraries like [FFTW](http://www.fftw.org/) and [libsndfile](http://www.mega-nerd.com/libsndfile/).

- Special thanks to the open-source community for tools like [OpenBLAS](https://www.openblas.net/), [libpng](http://www.libpng.org/pub/png/libpng.html), and [OpenMP](https://www.openmp.org/).

- Gratitude to [lucasb-eyer/libheatmap](https://github.com/lucasb-eyer/libheatmap) for the heatmap visualization module used in this project.

- Credit to [lieff/minimp3](https://github.com/lieff/minimp3) for the lightweight MP3 decoding library.

