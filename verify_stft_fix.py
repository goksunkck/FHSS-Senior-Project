import numpy as np
import sys
import os

# Add src to path
sys.path.append(os.path.join(os.getcwd(), 'src', 'python'))

from dsp_utils import compute_stft_matrix

# Test STFT scaling with a simple tone
fs = 10e6
duration = 0.001  # 1 ms
t = np.arange(0, duration, 1/fs)

# Create a simple 1 MHz tone with amplitude = 1
freq = 1e6
signal = np.exp(1j * 2 * np.pi * freq * t)

print(f"Signal length: {len(signal)}")
print(f"Signal power: {np.mean(np.abs(signal)**2):.6f}")
print(f"Signal power (dB): {10*np.log10(np.mean(np.abs(signal)**2)):.2f}")

# STFT parameters (matching our dataset)
window_len = 128
overlap_len = 96
nfft = 1024

# Use our function
log_mag = compute_stft_matrix(signal, fs, nfft, window_len, overlap_len)

print(f"\nPython STFT Output (with fix):")
print(f"  Shape: {log_mag.shape}")
print(f"  Max magnitude (dB): {log_mag.max():.2f}")
print(f"  Mean magnitude (dB): {log_mag.mean():.2f}")

print("\n" + "="*60)
print("COMPARISON:")
print("="*60)
print(f"MATLAB spectrogram peak: 36.72 dB")
print(f"Python STFT peak (fixed): {log_mag.max():.2f} dB")
print(f"Difference: {log_mag.max() - 36.72:.2f} dB")
