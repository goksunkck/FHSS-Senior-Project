import h5py
import numpy as np

"""
Re-check: The user is right - within ONE signal, there should be 
a deterministic Gold code pattern. Let me verify this.
"""

metadata_file = 'data/synthetic/prepared_prediction_sequences_python.h5'
dataset_file = 'data/synthetic/classification_dataset_stft_random_deg10_python.h5'

with h5py.File(metadata_file, 'r') as f:
    seq_starts_train = f['sequence_starts_train'][:]
    labels_train = f['YTrain'][:]
    lookback = f.attrs['lookback_window']

with h5py.File(dataset_file, 'r') as f:
    Y_data = f['Y_data']
    
    # Get the FIRST signal (should be hops 0-255)
    signal_0_labels = Y_data[0:256, 0]
    
    print("="*60)
    print("Checking if Gold code is deterministic within ONE signal")
    print("="*60)
    print(f"\nSignal 0 (first 256 hops):")
    print(f"  Label sequence (first 30): {signal_0_labels[:30]}")
    
    # Check autocorrelation of THIS signal's labels
    from scipy.signal import correlate
    
    autocorr = correlate(signal_0_labels - signal_0_labels.mean(), 
                        signal_0_labels - signal_0_labels.mean(), 
                        mode='full')
    autocorr = autocorr[len(autocorr)//2:]  # Keep only positive lags
    autocorr = autocorr / autocorr[0]  # Normalize
    
    print(f"\nAutocorrelation of Signal 0's hop sequence:")
    print(f"  Lag 0: {autocorr[0]:.4f}")
    print(f"  Lag 1: {autocorr[1]:.4f}")
    print(f"  Lag 2: {autocorr[2]:.4f}")
    print(f"  Lag 10: {autocorr[10]:.4f}")
    
    # Now check: Given hops [t-35 : t], can we determine hop t+1?
    # Let's check if the pattern is unique
    print("\n" + "="*60)
    print("Checking if lookback=35 uniquely determines next hop")
    print("="*60)
    
    # Create all possible 35-hop sequences from this signal
    sequences = {}
    for i in range(256 - lookback):
        pattern = tuple(signal_0_labels[i : i + lookback])
        next_hop = signal_0_labels[i + lookback]
        
        if pattern in sequences:
            if sequences[pattern] != next_hop:
                print(f"\n⚠️  COLLISION! Same 35-hop pattern maps to different next hops:")
                print(f"  Pattern ending at hop {i}: next = {next_hop}")
                print(f"  Previously seen: next = {sequences[pattern]}")
        else:
            sequences[pattern] = next_hop
    
    print(f"\nUnique 35-hop patterns in signal 0: {len(sequences)}")
    print(f"Total sequences checked: {256 - lookback}")
    
    if len(sequences) == 256 - lookback:
        print("\n✓ Good! Each 35-hop pattern uniquely determines the next hop.")
        print("  Within this signal, the task IS learnable.")
    
    # But now check: Are these patterns shared across different signals?
    print("\n" + "="*60)
    print("Checking if patterns generalize across signals")
    print("="*60)
    
    # Get Signal 1
    signal_1_labels = Y_data[256:512, 0]
    
    print(f"\nSignal 1 (hops 256-511):")
    print(f"  Label sequence (first 30): {signal_1_labels[:30]}")
    
    # Check if ANY of signal 0's patterns appear in signal 1
    signal_1_sequences = set()
    for i in range(256 - lookback):
        pattern = tuple(signal_1_labels[i : i + lookback])
        signal_1_sequences.add(pattern)
    
    overlap = set(sequences.keys()) & signal_1_sequences
    print(f"\nPatterns in Signal 0: {len(sequences)}")
    print(f"Patterns in Signal 1: {len(signal_1_sequences)}")
    print(f"Overlap: {len(overlap)}")
    
    if len(overlap) == 0:
        print("\n⚠️  ZERO overlap! Each signal has completely unique patterns.")
        print("  The model CANNOT generalize across signals.")
        print("  It would need to memorize 1100 * 221 = 243K unique patterns!")

print("\n" + "="*60)
print("CONCLUSION")
print("="*60)
