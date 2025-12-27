import numpy as np
import scipy.signal

def fsk_modulate(bits, sim_params):
    """
    Simple Frequency Shift Keying Modulator.
    FSK: Map bits to frequencies.
    Input:
        bits: binary array (0/1)
        sim_params: dictionary with keys:
            fs, symbolRate, freqSeparation, M, samplesPerSymbol
    Output:
        modulated_signal: Complex baseband signal? Or Passband?
        MATLAB script generates 'modulatedData' which seems to be complex baseband 
        before being upconverted by the hopping carrier.
        Wait, MATLAB's comm.FSKModulator output is Baseband?
        Usually yes. 
        Freq Separation is delta_f. 
        M=2 (Binary FSK).
        Frequencies usually: -sep/2, +sep/2 (for 2-FSK centered).
        Or 0, sep? 
        MATLAB FSKModulator default: 
        "The object maps the input symbols to valid frequencies... 
        Frequency separation is the distance between adjacent frequencies."
        Usually deviations are +/- (M-1)*sep/2 ...
        
        Let's implement generating a phase ramp.
        freqs = (2*symbol - (M-1)) * (freqSeparation / 2)
        for M=2, symbols 0, 1:
        0 -> (0 - 1) * sep/2 = -sep/2
        1 -> (2 - 1) * sep/2 = +sep/2
    """
    M = sim_params['M']
    samples_per_symbol = sim_params['samplesPerSymbol']
    fs = sim_params['fs']
    freq_sep = sim_params['freqSeparation']
    
    # 1. Map bits to symbols (0..M-1)
    # Since M=2, bits are symbols.
    symbols = bits.copy() # (N, )
    
    # 2. Map symbols to frequencies
    # Formula for M-FSK frequencies relative to carrier:
    # f_i = (2*i - (M-1)) * (freq_sep / 2)
    # i in 0..M-1
    # For M=2: 0 -> -sep/2, 1 -> +sep/2
    freqs = (2 * symbols - (M - 1)) * (freq_sep / 2.0)
    
    # 3. Generate phase
    # We need to repeat each frequency for 'samples_per_symbol'
    # But phase must be continuous for CPM/CPFSK? 
    # MATLAB comm.FSKModulator defaults to "Continuous Phase = True".
    # So we integrate frequency to get phase.
    
    # Expand freqs to sample rate
    freqs_upsampled = np.repeat(freqs, samples_per_symbol)
    
    # Integrate: phase = 2*pi * cumsum(freq/fs)
    # delta_t = 1/fs
    phase = 2 * np.pi * np.cumsum(freqs_upsampled) / fs
    
    modulated = np.exp(1j * phase)
    return modulated

def add_awgn(signal, snr_db):
    """
    Add AWGN noise to signal for a given SNR (dB), based on measured signal power.
    """
    sig_power = np.mean(np.abs(signal)**2)
    sig_power_db = 10 * np.log10(sig_power)
    
    noise_db = sig_power_db - snr_db
    noise_power = 10 ** (noise_db / 10)
    
    # Complex noise
    # Power = E[|n|^2] = E[nr^2 + ni^2] = 2 * sigma^2
    # So sigma = sqrt(noise_power / 2)
    sigma = np.sqrt(noise_power / 2)
    
    noise = sigma * (np.random.randn(len(signal)) + 1j * np.random.randn(len(signal)))
    return signal + noise

def compute_stft_matrix(signal_chunk, fs, nfft, window_len, overlap_len):
    """
    Compute STFT matching MATLAB's spectrogram 'centered' output.
    Output: (Freq, Time) matrix, Log-Magnitude.
    """
    # scipy.signal.spectrogram returns (f, t, Sxx)
    # Sxx is usually Power or Magnitude squared if 'spectrum'? 
    # Let's use stft to get complex first equivalent to MATLAB's [S,F,T] = spectrogram(...)
    # MATLAB: spectrogram(x, window, overlap, nfft, fs, 'centered')
    
    win = scipy.signal.get_window('hamming', window_len)
    
    # return_onesided=False gives two-sided.
    # verify nperseg vs nfft.
    f, t, Zxx = scipy.signal.stft(
        signal_chunk, 
        fs=fs, 
        window=win, 
        nperseg=window_len, 
        noverlap=overlap_len, 
        nfft=nfft, 
        return_onesided=False,
        boundary=None, # MATLAB default often behavior at edges differs, 'centered' in MATLAB spectrogram refers to FREQUENCY axis centering? Yes.
        padded=False # MATLAB might truncate or pad differently.
    )
    
    # Python returns 0..fs. MATLAB 'centered' returns -fs/2..fs/2.
    # We need fftshift to center it.
    Zxx = np.fft.fftshift(Zxx, axes=0)
    
    # CRITICAL MATLAB COMPATIBILITY FIX:
    # MATLAB's spectrogram() applies a scaling factor that scipy.signal.stft does not.
    # The scaling factor is approximately: window_len * coherent_gain
    # For Hamming window (128 length): sum(hamming(128))/128 ≈ 0.54
    # So: 128 * 0.54 ≈ 69.12
    # This matches MATLAB's output: 20*log10(69) ≈ 36.78 dB
    win_array = scipy.signal.get_window('hamming', window_len)
    coherent_gain = np.sum(win_array) / len(win_array)
    scaling_factor = window_len * coherent_gain
    Zxx = Zxx * scaling_factor
    
    # Log-Mag
    mag = np.abs(Zxx)
    log_mag = 20 * np.log10(mag + np.finfo(float).eps)
    
    # IMPORTANT: Check shape.
    # MATLAB: (Freq=1024, Time=3). 
    # Python: (Freq=1024, Time=K).
    # We need to ensure K=3.
    # samplesPerHop is calculated in main. If lengths match, it should work.
    
    return log_mag.astype(np.float32)
