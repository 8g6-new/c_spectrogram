import librosa
import numpy as np
import cv2
import os

def normalize_data(data):
    return (data - np.min(data)) / (np.max(data) - np.min(data))

def compute_mel(file_path, op_path,sr=44100, n_fft=4096, hop_length=64, n_mels=512, fmin=0, fmax=10000):
    y, sr               = librosa.load(file_path, sr=sr, dtype=np.float32)
    mel_spec            = librosa.feature.melspectrogram(y=y, sr=sr, n_fft=n_fft, hop_length=hop_length, n_mels=n_mels, fmin=fmin, fmax=fmax)
    mel_spec_db         = librosa.amplitude_to_db(mel_spec, ref=np.max)
    mel_spec_normalized = normalize_data(mel_spec_db)
    mel_spec_flipped    = np.flip(mel_spec_normalized, axis=0)
    mel_colormap        = cv2.applyColorMap((mel_spec_flipped * 255).astype(np.uint8), cv2.COLORMAP_JET)

    
    cv2.imwrite(op_path, mel_colormap)

file_path = './tests/files/black_woodpecker.wav'
op_path   = './python_test.png'
compute_mel(file_path, op_path)
