# 🎧 CARA (C Acoustic Representation & Analysis): High-Performance Audio Signal Processing and Visualization Pipeline

**CARA** is a high-performance C library for audio signal processing and visualization, featuring Short-Time Fourier Transform (STFT), Mel spectrograms, Mel-Frequency Cepstral Coefficients (MFCC), and professional-grade heatmap visualizations. Optimized for large-scale audio datasets, it leverages [FFTW](http://www.fftw.org/) with wisdom caching, [OpenMP](https://www.openmp.org/) parallelization, and BLAS ([OpenBLAS](https://www.openblas.net/)) for fast matrix operations. The library supports multiple audio formats (WAV, FLAC, MP3) via [libsndfile](https://libsndfile.github.io/libsndfile/) and [minimp3](https://github.com/lieff/minimp3), and offers customizable visualizations with extensive color schemes.

## Key Features

- **Audio I/O**: Reads WAV, AAC, MP3 etc with auto-detection. [minimp3](https://github.com/lieff/minimp3) decodes MP3s as fast as minimp3 allows and for other formats it use [libsndfile](https://libsndfile.github.io/libsndfile/).
- **Short-Time Fourier Transform (STFT)**: Uses FFTW with wisdom caching to plan FFTs, but it’s still slower than Librosa for reasons I can’t quite pin down. Supports window functions like Hann, Hamming, Blackman, and more. Tweak window size, hop size, and frequency bounds to your heart’s content.
- **Mel Spectrogram**: Builds Mel filter banks dynamically and tries to go fast with BLAS (`cblas_sdot`) and OpenMP. it’s not fast enough, and the output went wonky when we pushed too hard. Still, it’s got optional dB scaling with branchless computation for a bit of cool.
- **Mel-Frequency Cepstral Coefficients (MFCC)**: Computes MFCCs with precomputed DCT coefficients and BLAS. OpenMP helps, but it’s not blowing anyone away. Heatmap visualizations look nice, though, with customizable colors.
- **Visualization**: Spits out STFT, Mel spectrograms, and MFCCs as PNG heatmaps via [libheatmap](https://github.com/lucasb-eyer/libheatmap). Tons of color schemes (Blues, Viridis, Inferno, Plasma) in discrete, mixed, mixed_exp, and soft flavors. Cache-friendly memory ops keep rendering snappy.
- **Benchmarking**: Microsecond-precision profiling for STFT, Mel, MFCC, and visualization. Ranked timing reports with color-coded bars show you exactly where the pipeline’s dragging its feet. Outputs JSON and raw formats for analysis.
- **Performance Optimizations**: OpenMP parallelizes everything, FFTW wisdom caching saves time, and BLAS tries to speed up matrix ops (with mixed results). Compiler flags (`-ffast-math`, `-march=native`, `-funroll-loops`, LTO) and aligned memory allocations push the hardware, not at Librosa’s level yet.
- **Applications**: Great for bioacoustics (e.g., bird calls with `tests/files/black_woodpecker.wav`, `tests/files/173.mp3`), large-scale audio processing, ML feature extraction, and DSP research.

## 💡 Motivation

The main motivation behind this project was to gain a deeper understanding of both **C** and **digital signal processing (DSP)**. While there are countless tutorials on how to **use** MFCCs and Mel filter banks, very few actually explain how to **compute** them from scratch. The process was often fragmented or hidden behind library calls.

When searching for minimalist MFCC pipelines, I came across excellent projects like [rust-mfcc](https://github.com/bytesnake/mfcc), which performed impressively — about **2.5× faster than Librosa** on synthetic benchmarks ([Colab Notebook](https://github.com/8g6-new/mfcc_rust_bench/blob/master/rust_vs_python.ipynb)).  
However, they often rely on external dependencies and abstractions that obscure what's happening under the hood.

I noticed a lack of **simple, dependency-free, well-structured C implementations** of STFT, Mel spectrograms, and MFCCs that emphasize:

1. **Readability** – Code that beginners in C can actually follow  
2. **Educational Value** – A step-by-step DSP pipeline laid bare  
3. **Transparency** – Each transform is explicitly written (FFT, Mel bank, DCT)

As I built this project, I came to understand and appreciate:
- How **windowing**, **hop size**, and **FFT resolution** interact  
- The inner workings of **Mel filter bank construction**  
- How to derive **MFCCs using DCT**, and why the coefficients matter  
- The performance implications of **memory layout**, **cache locality**, and **contiguous memory access**  
- How small details like **loop nesting**, **BLAS vectorization**, and **data alignment** can drastically affect speed

This project isn't **currently** trying to beat Librosa or Rust DSP libraries in performance — **though future optimizations may close the gap**.

Instead, it's meant to be a **clear, hackable, and minimalist reference** for students, hobbyists, and anyone who wants to learn DSP by building it from the ground up.

If it helps others demystify the DSP pipeline or write their own from scratch, then it's done its job.

## Pipeline Overview

The following diagram illustrates the audio processing and visualization pipeline:

```mermaid
%%{init: {
  "theme": "base",
  "themeVariables": {
    "primaryColor": "#1e1e1e",          %% Node background (very dark)
    "primaryTextColor": "#ffffff",      %% Text color (white)
    "primaryBorderColor": "#ffaa00",    %% Orange border for visibility
    "clusterBkg": "#2a2a2a",            %% Subgraph background (dark gray)
    "clusterBorder": "#ffaa00",         %% Subgraph border color
    "lineColor": "#ffaa00",             %% Arrow color
    "fontSize": "14px",
    "fontFamily": "monospace"
  }
}}%%

flowchart TD

    %% Input and Decoding
    A["📥 Audio Input (.wav / .mp3)"] --> B["🔍 Auto File Type Detection"]
    B --> C{"🧩 Format Type"}
    C -->|MP3| D["🎧 Decode with minimp3"]
    C -->|Other| E["🎵 Read with libsndfile"]
    D --> F["🎚️ Normalize → Float32"]
    E --> F

    %% Feature Extraction
    subgraph Feature Extraction
        F --> G["🪟 Apply Window Function (e.g., Hann)"]
        G --> H["⚡ STFT (FFTW + Wisdom)"]
        H --> I["📊 Extract Magnitudes & Phases"]
        I --> J["🎚️ Apply Mel Filter Bank (BLAS)"]
        J --> K["🎯 Compute MFCC (DCT)"]
    end

    %% Visualization
    subgraph Visualization
        H --> V1["🖼️ STFT Heatmap"]
        J --> V2["🎨 Mel Spectrogram"]
        K --> V3["🌡️ MFCC Heatmap"]
    end

    %% Benchmarking
    subgraph Benchmarking
        H --> B1["⏱️ Time STFT"]
        J --> B2["⏱️ Time Mel Computation"]
        K --> B3["⏱️ Time MFCC Extraction"]
        V1 --> B4["⏱️ Time Plot Generation"]
    end

```

## Performance Highlights

- **MP3 Decoding & PNG Saving**: Shockingly fast — faster than Librosa, faster than anything in Python. minimp3 and libpng just show up, do their job, and leave. If the whole pipeline moved like this, we’d be done before the coffee brewed.

- **STFT & Mel Spectrogram**: Still slower than Librosa — even with FFTW wisdom caching and OpenMP. Not sure why. Librosa somehow still beats it. The Mel spectrogram part was especially disappointing: I tried to make it fast with BLAS, but the output came out wrong. Only one loop could be vectorized — the other two just sat there, immune to optimization. The filter bank creation is clean, but the actual dot-product part still suffers under that cursed 2-level nested loop ( I could eliminate 1 loop via BALS though, kinda win ig).

- **Scalability**: OpenMP does help  not so much , very much unsable until you compare the core DSP to librosa

## Requirements

- **Compiler**: GCC or Clang with C11 support.
- **Dependencies**:
  - **FFTW3** ([FFTW](http://www.fftw.org/)) for fast Fourier transforms.
  - **libsndfile** ([libsndfile](https://libsndfile.github.io/libsndfile/)) for WAV/FLAC file handling.
  - **OpenMP** ([OpenMP](https://www.openmp.org/)) for parallel processing.
  - **BLAS** (e.g., [OpenBLAS](https://www.openblas.net/)) for matrix operations.
  - **libpng** ([libpng](http://www.libpng.org/pub/png/libpng.html)) for PNG output.
- **Hardware**: Modern CPU recommended for optimal performance.

## Installation

### Step 1: Install Dependencies
For Ubuntu/Debian:
```bash
sudo apt-get update
sudo apt-get install libfftw3-dev libsndfile1-dev libopenblas-dev libpng-dev libomp-dev
```
For OpenCV-like color schemes (optional):
```bash
sudo apt-get install libopencv-dev
```

### Step 2: Clone the Repository
```bash
git clone https://github.com/8g6-new/c_spectrogram
cd c_spectrogram
```

### Step 3: Build the Project
Choose a build target:
- **Built-in color schemes**:
  ```bash
  make builtin
  ```
- **OpenCV-like color schemes**:
  ```bash
  make opencv_like
  ```
- **Shared library**:
  ```bash
  make shared
  ```
- **Debug builds**:
  ```bash
  make debug_builtin
  make debug_opencv_like
  ```

The build creates executables in `build/builtin` or `build/opencv` and generates FFTW wisdom files in `cache/FFT` (e.g., `1024.wisdom`).

## Usage

### Command-Line Interface
Run the `main` program to process an audio file and generate STFT, Mel spectrogram, and MFCC visualizations:
```bash
./build/builtin/main <input_file> <output_prefix> <window_size> <hop_size> <window_type> <num_mel_banks> <min_mel> <max_mel> <num_mfcc_coeffs> <cs_stft> <cs_mel> <cs_mfcc> <cache_folder>
```

**Parameters**:
- `input_file`: Path to audio file (e.g., `tests/files/black_woodpecker.wav`).
- `output_prefix`: Prefix for output PNG files (e.g., `outputs/black_woodpecker`).
- `window_size`: STFT window size (e.g., 2048).
- `hop_size`: Hop size for STFT (e.g., 512).
- `window_type`: Window function (e.g., `hann`, `hamming`, `blackman`).
- `num_mel_banks`: Number of Mel filters (e.g., 40).
- `min_mel`, `max_mel`: Frequency range for Mel filters (e.g., 20.0, 8000.0).
- `num_mfcc_coeffs`: Number of MFCC coefficients (e.g., 13).
- `cs_stft`, `cs_mel`, `cs_mfcc`: Color scheme indices (e.g., 0 for default, see `output/colors.json`).
- `cache_folder`: Directory for FFTW wisdom files (e.g., `cache/FFT`).

**Example**:
```bash
./build/builtin/main tests/files/black_woodpecker.wav outputs/black_woodpecker 2048 512 hann 40 20.0 8000.0 13 0 0 0 cache/FFT
```

**Output**:
- PNG files: `outputs/black_woodpecker_stft.png`, `outputs/black_woodpecker_mel.png`, `outputs/black_woodpecker_mfcc.png`.
- Console output: Audio metadata (duration, sample rate) and ranked benchmark timings.

### Programmatic Usage
Below is a simplified example of using the library in C:

```c
#include "headers/audio_tools/audio_io.h"
#include "headers/audio_tools/spectral_features.h"
#include "headers/audio_tools/audio_visualizer.h"
#include "headers/utils/bench.h"

int main() {
    benchmark_init(); // Initialize benchmarking

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

    // Compute MFCC
    size_t num_coeffs = 13;
    mffc_t dft_coeffs = precompute_cosine_coeffs(num_filters, num_coeffs);

    // Set visualization bounds
    bounds2d_t bounds = { .freq = {20.0f, 8000.0f} };
    init_bounds(&bounds, &result);
    set_limits(&bounds, result.num_frequencies, result.output_size);

    // Copy STFT magnitudes
    float *contiguous_mem = malloc((bounds.freq.end_d - bounds.freq.start_d) * (bounds.time.end_d - bounds.time.start_d) * sizeof(float));
    fast_copy(contiguous_mem, result.magnitudes, &bounds, result.num_frequencies);

    // Visualize STFT spectrogram
    plot_t settings = { .db = true, .cs_enum = CS_Blues, .output_file = "outputs/stft.png" };
    settings.bg_color[0] = 0; settings.bg_color[1] = 0; settings.bg_color[2] = 0; settings.bg_color[3] = 255;
    spectrogram(contiguous_mem, &bounds, &settings);

    // Compute and visualize Mel spectrogram
    settings.output_file = "outputs/mel.png";
    float *mel_values = mel_spectrogram(contiguous_mem, num_filters, result.num_frequencies, mel_filter_bank, &bounds, &settings);

    // Compute and visualize MFCC
    settings.output_file = "outputs/mfcc.png";
    mfcc(mel_values, &dft_coeffs, &bounds, &settings);

    // Print benchmark results
    print_bench_ranked();

    // Cleanup
    free(mel_values);
    free(mel_filter_bank);
    free(contiguous_mem);
    free(window_values);
    free_stft(&result);
    free_fft(&fft);
    free_audio(&audio);
    free(dft_coeffs.coeffs);
    free(melbank.freq_indexs);
    free(melbank.weights);
    return 0;
}
```


## 📊 Full Function Visualizations

Below are visualizations produced by the pipeline for a single audio input (`black_woodpecker.wav`), using various color schemes and stages of the DSP pipeline: STFT, Mel Spectrogram, and MFCC.

---

### 🎛️ Function Outputs

| Output Type               | Color Scheme     | Description                                         | Preview |
|---------------------------|------------------|-----------------------------------------------------|---------|
| **STFT Spectrogram**      | Cividis          | Raw STFT data visualized as a spectrogram           | ![STFT Cividis](outputs/functions/black_woodpecke_stft.png) |
| **Mel Spectrogram**       | Cividis          | Mel filter bank output from STFT magnitudes         | ![Mel Spectrogram](outputs/functions/black_woodpecke_mel.png) |
| **MFCC**                  | Blues Soft       | Cepstral coefficients                               | ![MFCC](outputs/functions/black_woodpecke_mfcc.png) |

---

### 🐢 STFT Spectrograms (Built-in Color Schemes)

| Colormap        | Style      | Preview |
|------------------|------------|---------|
| **Blues**        | Soft       | ![STFT Blues Soft](outputs/colorschemes/libheatmap_defaults/soft/black_woodpecker_stft_Blues_soft.png) |
| **Spectral**     | Discrete   | ![STFT Spectral Discrete](outputs/colorschemes/libheatmap_defaults/discrete/black_woodpecker_stft_Spectral_discrete.png) |

---

### 🎨 STFT Spectrograms (OpenCV-like Color Schemes)

| Colormap        | Description                                       | Preview |
|------------------|---------------------------------------------------|---------|
| **Viridis**      | Scientific colormap emphasizing smooth gradients | ![STFT Viridis](outputs/colorschemes/opencv_like/images/black_woodpecker_stft_Viridis.png) |
| **Jet**          | High-contrast legacy colormap                    | ![STFT Jet](outputs/colorschemes/opencv_like/images/black_woodpecker_stft_Jet.png) |


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
- `colorschemes/libheatmap_defaults`:
  - `discrete`: High-contrast colormaps (e.g., `black_woodpecker_stft_Blues_discrete.png`).
  - `mixed`: Smooth color transitions (e.g., `black_woodpecker_stft_Blues_mixed.png`).
  - `mixed_exp`: Exponentially scaled colors (e.g., `black_woodpecker_stft_Blues_mixed_exp.png`).
  - `soft`: Softened gradients (e.g., `black_woodpecker_stft_Blues_soft.png`).
- `colorschemes/opencv_like/images`: OpenCV-inspired colormaps (e.g., `black_woodpecker_stft_Viridis.png`).
- `functions`: Mel spectrograms and MFCC heatmaps (e.g., `black_woodpecker_mel.png`, `black_woodpecker_mfcc.png`).

## Benchmarking Output
The `print_bench_ranked` function generates a ranked table of execution times:
- Columns: Function name, execution time (µs, ms, or s), percentage of total runtime.
- Visual: Color-coded bars for quick bottleneck identification.

**Example**:
```
### 🔍 Sample Benchmark Output

```text
Running builtin...
Input Filename       : ./tests/files/173.mp3
Output Filename      : bird
Window Size          : 2048
Hop Size             : 128
Window Type          : hann
Number of Filters    : 128
Min Mel Frequency    : 0.00
Max Mel Frequency    : 7500.00
Number of Coeffs     : 64
Cache STFT Channels  : 2
Cache Mel Channels   : 3
Cache MFCC Channels  : 4
Cache Folder         : ./cache/FFT
./tests/files/173.mp3 auto detected to be audio/mpeg
duration:58.032
channels:1
num_samples:2785536
sample_rate:48000.0000
file_size_bytes:784356
Cache not found or import failed. Creating FFT plan...
Error saving wisdom to file: ./cache/FFT/2048.wisdom
enum 2 => libheatmap builtin color scheme : main type PRGn , subtype mixed
stft ended
Time bounds:
  Start (f): 0.00
  End   (f): 0.00
  Start (d): 0
  End   (d): 21747
Frequency bounds:
  Start (f): 0.00
  End   (f): 7500.00
  Start (d): 0
  End   (d): 320

copy ended

bird_stft.png saved
bird_mel.png saved
bird_mfcc.png saved
---------------------------------------------------------
| Function             | Exec Time    | % of total runtime |
---------------------------------------------------------
| stft_plot            |      7.206 s | 49.4490% |
[▰▰▰▰▰▰▰▰▰           ]
| mel                  |      3.518 s | 24.1430% |
[▰▰▰▰                ]
| mfcc                 |      2.283 s | 15.6638% |
[▰▰▰                 ]
| stft                 |   925.131 ms | 6.3480% |
[▰                   ]
| dec_mp3              |   482.269 ms | 3.3092% |
[                    ]
| fetch fft cache 1    |   129.612 ms | 0.8894% |
[                    ]
| copy                 |    27.156 ms | 0.1863% |
[                    ]
| file_read            |      658 µs | 0.0045% |
[                    ]
| dtft coff            |      391 µs | 0.0027% |
[                    ]
| mel filter bank      |      390 µs | 0.0027% |
[                    ]
| stft window          |      146 µs | 0.0010% |
[                    ]
| auto_det             |       65 µs | 0.0004% |
[                    ]
---------------------------------------------------------
```

Mel looks significantly slower because it calls additional weighted points from the libheatmap lib, which adds delay, doing 2 separate loops was found to be even slower
> **Note** : bechmarked in colab


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

- **Multi-threading**: Set `OMP_NUM_THREADS` to match CPU cores (e.g., `export OMP_NUM_THREADS=8`).
- **BLAS**: Use optimized BLAS (e.g., [OpenBLAS](https://www.openblas.net/)) for faster Mel and MFCC computations.
- **Wisdom Caching**: Pre-generate FFTW wisdom files for common window sizes.
- **Matrix Operations**: Replace `cblas_sdot` with `cblas_sgemm` for batched matrix multiplications in `mel_spectrogram` and `mfcc`.
- **Fast Math**: Test `-ffast-math` for numerical stability, as it may introduce inaccuracies.
- **GPU**: Pipeline is GPU-ready for STFT ([cuFFT](https://developer.nvidia.com/cufft)), Mel spectrograms ([cuBLAS](https://developer.nvidia.com/cublas)), and visualization (CUDA kernels).

## Future Work

- **Explicit SIMD Support**: Implement explicit SIMD optimizations (e.g., SSE, SSE2, AVX, AVX2) for STFT, Mel spectrogram, and MFCC computations, beyond current implicit support via `minimp3` and compiler flags.
- **GPU Acceleration**: Implement CUDA-based STFT ([cuFFT](https://developer.nvidia.com/cufft)), Mel spectrograms ([cuBLAS](https://developer.nvidia.com/cublas)), and visualization for Librosa-like performance.
- **Sparse Matrix Operations**: Use CSR format for Mel filter banks to reduce memory and computation.
- **Real-Time Processing**: Support streaming audio analysis.
- **Enhanced Benchmarking**: Add memory usage and CPU/GPU utilization metrics to `bench.h`.
- **Simple Heatmap**: Replace `heatmap_add_weighted_point()` in `libheatmap` with a custom heatmap generator for faster rendering.
- **Memory Pooling**: Implement memory pooling using a memory arena to improve utilization and prevent issues like use-after-free.
- **Advanced Memory Management**: Integrate buddy allocators or other techniques to optimize memory allocation and reduce fragmentation.
- **Documentation**: Add detailed API docs and usage examples to `headers/` and `README.md`.

## 📄 License

- **Code**: Licensed under the [MIT License](./LICENSE).  
  You are free to use, modify, and distribute the code, including for commercial purposes, with proper attribution.

Credit this work if used in research or applications.

## Acknowledgments

- Inspired by [Librosa](https://librosa.org/) for high-performance audio processing in C.
- Tested on bioacoustics datasets (e.g., bird calls), with thanks to [FFTW](http://www.fftw.org/) and [libsndfile](https://libsndfile.github.io/libsndfile/).
- Gratitude to [OpenBLAS](https://www.openblas.net/), [libpng](http://www.libpng.org/pub/png/libpng.html), and [OpenMP](https://www.openmp.org/).
- Thanks to [lucasb-eyer/libheatmap](https://github.com/lucasb-eyer/libheatmap) for the heatmap visualization module.
- Credit to [lieff/minimp3](https://github.com/lieff/minimp3) for lightweight MP3 decoding.
