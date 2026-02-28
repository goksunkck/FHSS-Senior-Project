
import numpy as np
import h5py
import matplotlib.pyplot as plt
import os
import sys

def main():
    print("Generating SNR Comparison Figure...")
    
    # We need to find 3 samples from the dataset: one at -10dB, one near 0dB, one near 10dB.
    # Dataset path
    h5_path = os.path.join('data', 'synthetic', 'm_sequence_dataset.h5') # Or gold
    
    if not os.path.exists(h5_path):
        print("Dataset not found!")
        return
        
    with h5py.File(h5_path, 'r') as f:
        X = f['X_data']
        SNR = f['SNR'][:]
        
        # Find indices
        # We want SNR close to -10, 0, 10
        # Since we generated -10, -8... they should be exact matches.
        
        idx_neg10 = np.where(np.isclose(SNR, -10.0))[0]
        idx_0 = np.where(np.isclose(SNR, 0.0))[0]
        idx_10 = np.where(np.isclose(SNR, 10.0))[0]
        
        if len(idx_neg10) == 0 or len(idx_0) == 0 or len(idx_10) == 0:
            print("Could not find required SNRs.")
            return
            
        # Pick random sample from each
        s_neg10 = X[idx_neg10[0]]
        s_0 = X[idx_0[0]]
        s_10 = X[idx_10[0]]
        
        # Plot
        fig, axes = plt.subplots(1, 3, figsize=(15, 5))
        
        # Min/Max for unified color scale (optional, but better to normalize per plot for visibility)
        # Actually, let's just plot them individually.
        
        # Transpose/Flip for proper display? 
        # Shape (1024, 3). Freq x Time.
        # Imshow expects (Height, Width).
        # Origin lower.
        
        im0 = axes[0].imshow(s_neg10, aspect='auto', origin='lower', cmap='jet')
        axes[0].set_title('SNR -10 dB')
        axes[0].set_xlabel('Time Hop')
        axes[0].set_ylabel('Frequency Bin')
        
        im1 = axes[1].imshow(s_0, aspect='auto', origin='lower', cmap='jet')
        axes[1].set_title('SNR 0 dB')
        axes[1].set_xlabel('Time Hop')
        
        im2 = axes[2].imshow(s_10, aspect='auto', origin='lower', cmap='jet')
        axes[2].set_title('SNR 10 dB')
        axes[2].set_xlabel('Time Hop')
        
        plt.tight_layout()
        
        out_path = 'results/snr_spectrograms.png'
        if not os.path.exists('results'): os.makedirs('results')
        plt.savefig(out_path, dpi=300)
        print(f"Saved {out_path}")

if __name__ == "__main__":
    main()
