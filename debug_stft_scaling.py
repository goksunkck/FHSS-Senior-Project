import numpy as np
import scipy.signal
import matplotlib.pyplot as plt

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

win = scipy.signal.get_window('hamming', window_len)

# Compute STFT
f, t_stft, Zxx = scipy.signal.stft(
    signal,
    fs=fs,
    window=win,
    nperseg=window_len,
    noverlap=overlap_len,
    nfft=nfft,
    return_onesided=False,
    boundary=None,
    padded=False
)

# Shift frequencies
Zxx = np.fft.fftshift(Zxx, axes=0)
f = np.fft.fftshift(f)

# Log magnitude
mag = np.abs(Zxx)
log_mag = 20 * np.log10(mag + np.finfo(float).eps)

print(f"\nSTFT Output:")
print(f"  Shape: {Zxx.shape}")
print(f"  Max magnitude (linear): {mag.max():.6f}")
print(f"  Max magnitude (dB): {log_mag.max():.2f}")
print(f"  Mean magnitude (dB): {log_mag.mean():.2f}")
print(f"  Mag at peak bin (dB): {log_mag[np.unravel_index(log_mag.argmax(), log_mag.shape)]:.2f}")

# Check window power
window_power = np.sum(win**2) / len(win)
print(f"\nWindow:")
print(f"  Power (normalized): {window_power:.6f}")
print(f"  Power (dB): {10*np.log10(window_power):.2f}")

# Expected: For a unit-power signal through a windowed STFT,
# the peak bin should be roughly proportional to:
# Amplitude * WindowLength * sqrt(WindowPower)
# For Hamming window, coherent gain ~ sum(window)/N ~ 0.54
coherent_gain = np.sum(win) / len(win)
print(f"  Coherent gain: {coherent_gain:.6f}")
print(f"  Coherent gain (dB): {20*np.log10(coherent_gain):.2f}")

print("\n" + "="*60)
print("EXPECTED vs ACTUAL:")
print("="*60)
print(f"Input signal power: 0 dB (unit amplitude complex tone)")
print(f"STFT peak magnitude: {log_mag.max():.2f} dB")
print(f"Difference: {log_mag.max() - 0:.2f} dB")
print("\nThis difference is the STFT processing gain.")
