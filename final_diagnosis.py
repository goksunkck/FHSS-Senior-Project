"""
Final diagnosis: Check if the STFT frames themselves encode
enough information to predict the next hop.

The issue: Each STFT frame just shows "signal is at frequency X".
But to predict the NEXT frequency, you need to know:
1. The current state of the Gold code LFSR
2. Which bit transition happens next

The STFT doesn't encode this information - it's just frequency content!
"""

import h5py
import numpy as np

metadata_file = 'data/synthetic/prepared_prediction_sequences_python.h5'
dataset_file = 'data/synthetic/classification_dataset_stft_random_deg10_python.h5'

with h5py.File(metadata_file, 'r') as f:
    seq_starts_train = f['sequence_starts_train'][:]
    labels_train = f['YTrain'][:]
    lookback = f.attrs['lookback_window']

with h5py.File(dataset_file, 'r') as f:
    Xdata = f['X_data']
    Ydata = f['Y_data']
    
    # Get one full signal (256 hops)
    signal_start = 0
    signal_hops = 256
    
    full_signal_X = Xdata[signal_start : signal_start + signal_hops]  # (256, 1024, 3)
    full_signal_Y = Ydata[signal_start : signal_start + signal_hops]  # (256,)
    
    print("="*60)
    print("CRITICAL QUESTION: Can we predict hop N+1 from hop N?")
    print("="*60)
    
    # For each hop, find the peak frequency bin
    peak_bins = []
    for i in range(len(full_signal_Y)):
        frame = full_signal_X[i]  # (1024, 3)
        # Average over time dimension
        freq_profile = frame.mean(axis=1)
        peak_bin = np.argmax(freq_profile)
        peak_bins.append(peak_bin)
    
    # Check: Does peak bin correlate with the label?
    # Label should tell us which of 8 frequencies
    print(f"\nFirst 20 hops:")
    print(f"  Labels (frequency index 0-7): {full_signal_Y[:20, 0]}")
    print(f"  Peak bins (out of 1024): {peak_bins[:20]}")
    
    # Check autocorrelation of labels (the Gold code sequence)
    from scipy.signal import correlate
    
    label_sequence = full_signal_Y[:, 0]
    autocorr = correlate(label_sequence - label_sequence.mean(), 
                        label_sequence - label_sequence.mean(), 
                        mode='full')
    autocorr = autocorr[len(autocorr)//2:]  # Keep only positive lags
    autocorr = autocorr / autocorr[0]  # Normalize
    
    print(f"\nAutocorrelation of Gold code sequence (labels):")
    print(f"  Lag 0: {autocorr[0]:.4f}")
    print(f"  Lag 1: {autocorr[1]:.4f}")
    print(f"  Lag 2: {autocorr[2]:.4f}")
    print(f"  Lag 5: {autocorr[5]:.4f}")
    print(f"  Lag 10: {autocorr[10]:.4f}")
    
    if abs(autocorr[1]) < 0.1:
        print("\n⚠️  CRITICAL ISSUE: Gold code is nearly white noise!")
        print("   Lag-1 autocorrelation ~0 means hop N gives NO info about hop N+1.")
        print("   This is expected for a good spreading code, but makes prediction impossible")
        print("   from frequency content alone!")
    
    print("\n" + "="*60)
    print("DIAGNOSIS:")
    print("="*60)
    print("The task as designed is IMPOSSIBLE without additional information.")
    print("")
    print("Why: Gold codes are designed to be pseudorandom. Each hop frequency")
    print("     depends on the LFSR state, but the STFT only shows 'signal at freq X'.")
    print("     There's no way to infer the LFSR state from frequency content.")
    print("")
    print("Solutions:")
    print("  1. Include the LFSR state or hop index as an auxiliary input")
    print("  2. Use a FIXED Gold code (same seed every time) so the model can memorize")
    print("  3. Change task: predict frequency from PREVIOUS hops (regression)")
    print("  4. Use phase information or chip-level timing (very complex)")
    print("="*60)
