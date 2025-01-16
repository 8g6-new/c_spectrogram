# c_spectrogram

<a href="https://wakatime.com/badge/user/018c3e96-0fae-432d-add6-28d53961f8b4/project/54fb69e3-a9fa-4548-b5ef-6494425e2ffb"><img src="https://wakatime.com/badge/user/018c3e96-0fae-432d-add6-28d53961f8b4/project/54fb69e3-a9fa-4548-b5ef-6494425e2ffb.svg" alt="wakatime"></a>


This repository contains tools and functions for audio signal processing and visualization, including functionality for performing Short-Time Fourier Transform (STFT), generating Mel spectrograms, applying different window functions, and rendering heatmaps.

## Features

- **STFT**: Compute the Short-Time Fourier Transform of an audio signal.
- **Window Functions**: Apply various windowing functions such as Hann, Hamming, Blackman, Bartlett, etc.
- **Mel Spectrogram**: Compute and visualize Mel spectrograms using Mel filter banks.
- **Heatmap Visualization**: Generate and save spectrograms as heatmaps in PNG format.

## Dependencies

- [FFTW](http://www.fftw.org/) - A library for performing fast Fourier transforms.
- [libsndfile](http://www.mega-nerd.com/libsndfile/) - A library for reading and writing sound files.
- [Heatmap](https://github.com/user/heatmap) - A custom heatmap rendering library (adjust link to your actual heatmap library).
- [PNG Tools](https://github.com/user/png_tools) - Custom PNG generation library (adjust link to your actual PNG tools library).

## Installation

1. Clone this repository:

   ```bash
   git clone https://github.com/yourusername/audio-signal-processing.git
   cd audio-signal-processing

