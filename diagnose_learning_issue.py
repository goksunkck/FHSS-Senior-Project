import h5py
import numpy as np
import matplotlib.pyplot as plt
import torch

# Load a few sequences and check if patterns are visible
metadata_file = 'data/synthetic/prepared_prediction_sequences_python.h5'
dataset_file = 'data/synthetic/classification_dataset_stft_random_deg10_python.h5'

with h5py.File(metadata_file, 'r') as f:
    seq_starts_train = f['sequence_starts_train'][:]
    labels_train = f['YTrain'][:]  # Not mapped
    lookback = f.attrs['lookback_window']
    
    # Get dataset filename
    ds_file = f.attrs['dataset_filename']
    
print(f"Dataset: {ds_file}")
print(f"Lookback window: {lookback}")
print(f"Train sequences: {len(seq_starts_train)}")
print(f"Unique labels: {np.unique(labels_train)}")
print(f"Label distribution: {np.bincount(labels_train)}")

# Load raw data
with h5py.File(dataset_file, 'r') as f:
    X_data = f['X_data']
    Y_data = f['Y_data']
    
    # Get norm stats
    if 'global_min' in f.attrs:
        g_min = float(f.attrs['global_min'])
        g_max = float(f.attrs['global_max'])
    else:
        g_min = float(f['global_min'][0])
        g_max = float(f['global_max'][0])
    
    print(f"\nNormalization: Min={g_min:.2f}, Max={g_max:.2f}")
    g_range = g_max - g_min
    
    # Sample 5 sequences, one from each class (if possible)
    unique_labels = np.unique(labels_train)
    print(f"\nSampling one sequence per class...")
    
    fig, axes = plt.subplots(len(unique_labels), 3, figsize=(15, 3*len(unique_labels)))
    
    for class_idx, target_class in enumerate(unique_labels):
        # Find a sequence with this label
        matching_indices = np.where(labels_train == target_class)[0]
        if len(matching_indices) == 0:
            continue
        
        sample_idx = matching_indices[0]
        start_idx = seq_starts_train[sample_idx]
        
        # Load sequence (lookback frames)
        sequence = X_data[start_idx : start_idx + lookback]
        
        # Normalize
        sequence = (sequence - g_min) / g_range
        
        # Plot first, middle, last frame
        frames_to_plot = [0, lookback//2, lookback-1]
        
        for plot_idx, frame_idx in enumerate(frames_to_plot):
            ax = axes[class_idx, plot_idx]
            
            # Frame shape: (1024, 3)
            # Show first time slice (frequency vs time)
            img = sequence[frame_idx]  # (1024, 3)
            
            ax.imshow(img, aspect='auto', cmap='viridis', vmin=0, vmax=1)
            ax.set_title(f"Class {target_class}, Frame {frame_idx}/{lookback}")
            ax.set_ylabel("Frequency Bin")
            ax.set_xlabel("Time")
            
            # Print stats
            print(f"  Class {target_class}, Frame {frame_idx}: "
                  f"Min={img.min():.3f}, Max={img.max():.3f}, Mean={img.mean():.3f}")
    
    plt.tight_layout()
    plt.savefig('results/diagnostic_samples_per_class.png', dpi=150)
    print("\nSaved: results/diagnostic_samples_per_class.png")
    
    # Now check: For a SINGLE sequence, do consecutive frames look similar?
    # (They should, because it's the same Gold code pattern)
    print("\n" + "="*60)
    print("Checking temporal consistency...")
    print("="*60)
    
    sample_idx = 0
    start_idx = seq_starts_train[sample_idx]
    label = labels_train[sample_idx]
    
    sequence = X_data[start_idx : start_idx + lookback]
    sequence = (sequence - g_min) / g_range
    
    # Compute correlation between consecutive frames
    correlations = []
    for i in range(len(sequence) - 1):
        corr = np.corrcoef(sequence[i].flatten(), sequence[i+1].flatten())[0, 1]
        correlations.append(corr)
    
    print(f"Correlation between consecutive frames:")
    print(f"  Mean: {np.mean(correlations):.4f}")
    print(f"  Std: {np.std(correlations):.4f}")
    print(f"  Min: {np.min(correlations):.4f}")
    print(f"  Max: {np.max(correlations):.4f}")
    
    if np.mean(correlations) < 0.1:
        print("\n⚠️  WARNING: Frames are nearly uncorrelated!")
        print("   This suggests each hop is independent, making temporal learning impossible.")
    
    # Check: Are different classes actually different?
    print("\n" + "="*60)
    print("Checking inter-class differences...")
    print("="*60)
    
    # Get one sequence from each class
    class_samples = []
    for target_class in unique_labels[:3]:  # Just check first 3
        matching_indices = np.where(labels_train == target_class)[0]
        if len(matching_indices) > 0:
            start_idx = seq_starts_train[matching_indices[0]]
            seq = X_data[start_idx : start_idx + lookback]
            seq = (seq - g_min) / g_range
            class_samples.append((target_class, seq))
    
    # Compute pairwise distances
    for i in range(len(class_samples)):
        for j in range(i+1, len(class_samples)):
            class_i, seq_i = class_samples[i]
            class_j, seq_j = class_samples[j]
            
            # Mean absolute difference
            diff = np.mean(np.abs(seq_i - seq_j))
            print(f"Class {class_i} vs Class {class_j}: Mean diff = {diff:.4f}")
    
    if diff < 0.01:
        print("\n⚠️  WARNING: Classes look nearly identical!")
        print("   The model has no pattern to learn.")

print("\n" + "="*60)
print("Diagnosis complete. Check results/diagnostic_samples_per_class.png")
print("="*60)
