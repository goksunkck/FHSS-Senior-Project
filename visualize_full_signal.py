import numpy as np
import h5py
import matplotlib.pyplot as plt
import os
import random

def visualize_full():
    print("--- Visualize Full Signal ---")
    
    dataset_path = 'data/synthetic/classification_dataset_stft_random_deg10_python.h5'
    
    if not os.path.exists(dataset_path):
        print(f"Error: Dataset {dataset_path} not found.")
        return

    with h5py.File(dataset_path, 'r') as f:
        X_data = f['X_data'] # Shape (TotalHops, 1024, 3)
        Y_data = f['Y_data']
        num_hops = int(f.attrs['numHopsPerSignal'])
        
        total_examples = X_data.shape[0]
        num_signals = total_examples // num_hops
        
        print(f"Dataset has {total_examples} hops, {num_signals} signals.")
        print(f"Hops per signal: {num_hops}")
        
        # Pick a random signal
        sig_idx = random.randint(0, num_signals - 1)
        print(f"Visualizing Signal Index: {sig_idx}")
        
        start_idx = sig_idx * num_hops
        end_idx = start_idx + num_hops
        
        # Read all hops for this signal
        # Shape: (256, 1024, 3)
        hops = X_data[start_idx:end_idx] 
        labels = Y_data[start_idx:end_idx].flatten()
        
        # Concatenate Time Dimension
        # We want to stack them horizontally.
        # Each hop is (1024, 3).
        # We want (1024, 3 * 256).
        
        # Transpose to (Freq, Time2)
        # hops[i] is (1024, 3).
        # We want to stack along axis 1.
        
        # Method 1: list of arrays
        hops_list = [hops[i] for i in range(num_hops)]
        full_spectrogram = np.hstack(hops_list) # (1024, 3*256)
        
        # Plot
        plt.figure(figsize=(20, 8))
        
        # Spectrogram
        # extent=[TimeStart, TimeEnd, FreqStart, FreqEnd]
        # Just use indices
        plt.imshow(full_spectrogram, aspect='auto', origin='lower', cmap='viridis')
        plt.title(f"Full Spectrogram - Signal {sig_idx} (SNR unknown)")
        plt.xlabel("Time (Hops * Slots)")
        plt.ylabel("Frequency Bins")
        plt.colorbar(label='Log Magnitude (dB)')
        
        # Overlay standard hopping pattern?
        # Labels are frequency indices (0..channels).
        # We need to map labels to Frequency Bins.
        # This mapping depends on grid size.
        # STFT 1024 bins. fs=10e6.
        # Hopset: 2e6 + k*spacing.
        # We can just plot the label indices on a secondary axis or overlay scatter.
        # Labels are 0..7. Frequency 0..7.
        # This doesn't map directly to 0..1024 bins without calculation.
        # But we can see the "stairs" pattern.
        
        # Let's save it.
        if not os.path.exists('results'):
            os.makedirs('results')
        
        save_path = 'results/full_signal_spectrogram.png'
        plt.tight_layout()
        plt.savefig(save_path)
        print(f"Saved full spectrogram to {save_path}")
        
        # Also plot the labels to verify pattern shape
        plt.figure(figsize=(20, 4))
        plt.plot(labels, 'r.-')
        plt.title(f"Hopping Pattern (Labels) - Signal {sig_idx}")
        plt.xlabel("Hop Index")
        plt.ylabel("Channel Index")
        plt.grid(True)
        plt.savefig('results/full_signal_labels.png')
        print(f"Saved labels plot to results/full_signal_labels.png")

if __name__ == "__main__":
    visualize_full()
