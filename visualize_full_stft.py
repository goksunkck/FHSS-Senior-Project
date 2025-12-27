
import numpy as np
import scipy.io
import h5py
import matplotlib.pyplot as plt
import os
import argparse

def load_mat_file(filepath, variable_name):
    """Loads a specific variable from a .mat file."""
    try:
        return scipy.io.loadmat(filepath)
    except NotImplementedError:
        import h5py
        with h5py.File(filepath, 'r') as f:
            data = {}
            for k, v in f.items():
                if isinstance(v, h5py.Dataset):
                    value = v[()]
                    if isinstance(value, np.ndarray):
                        if value.shape == (1, 1):
                            data[k] = value[0, 0]
                        else:
                            data[k] = value
                    else:
                        data[k] = value
            return data

def visualize_full_signal(target_snr, colormap='viridis'):
    print("--- 1. Configuration ---")
    # Parameters matching the generation script
    snr_levels_db = np.arange(-10, 12, 2) # -10, -8, ..., 10
    num_hopsets = 10
    num_signals_per_config = 3
    num_hops_per_signal = 256
    
    total_signals_per_snr = num_hopsets * num_signals_per_config
    total_hops_per_snr = total_signals_per_snr * num_hops_per_signal
    
    if target_snr not in snr_levels_db:
        print(f"Error: SNR {target_snr} dB not found in dataset.")
        print(f"Available SNRs: {snr_levels_db}")
        return

    # Calculate index offset
    snr_idx = np.where(snr_levels_db == target_snr)[0][0]
    
    # Start index in X_data (1-based in MATLAB, so 0-based here)
    # The hops are stored sequentially by SNR loop, then Hopset loop, then Signal loop.
    # We want the FIRST signal of the requested SNR.
    start_hop_idx = snr_idx * total_hops_per_snr
    end_hop_idx = start_hop_idx + num_hops_per_signal # Just one signal
    
    print(f"Target SNR: {target_snr} dB")
    print(f"Loading hops {start_hop_idx} to {end_hop_idx}...")

    print("--- 2. Load Dataset ---")
    dataset_path = 'data/synthetic/classification_dataset_stft_awgn.mat'
    
    if not os.path.exists(dataset_path):
        print(f"Error: Dataset not found at {dataset_path}")
        return

    try:
        data_contents = load_mat_file(dataset_path, 'X_data')
        X_data_in_memory = data_contents['X_data']
        print("Dataset loaded.")
    except Exception as e:
        print(f"Failed to load dataset: {e}")
        return

    # Handle HDF5 references if necessary
    if hasattr(X_data_in_memory, 'dtype') and X_data_in_memory.dtype == h5py.ref_dtype:
        print("Resolving HDF5 object references...")
        resolved_x_data = []
        with h5py.File(dataset_path, 'r') as f:
            for i in range(start_hop_idx, end_hop_idx):
                ref = X_data_in_memory.flat[i]
                resolved_x_data.append(f[ref][()])
        X_data_slice = resolved_x_data
    else:
        # Standard load
        X_data_slice = X_data_in_memory[start_hop_idx:end_hop_idx].flatten()

    print("--- 3. Construct Full Spectrogram ---")
    hop_frames = []
    
    for chunk in X_data_slice:
        if isinstance(chunk, np.ndarray) and chunk.dtype == object and chunk.size == 1:
            chunk = chunk[0]
            
        frame = chunk
        if frame.ndim == 3:
             frame = frame[:, :, 0] 
        
        if frame.shape[0] == 3 and frame.shape[1] > 3:
            frame = frame.T # Transpose to (Freq, Time)
        
        hop_frames.append(frame)

    full_spectrogram = np.concatenate(hop_frames, axis=1)
    print(f"Full spectrogram shape: {full_spectrogram.shape}")

    print("--- 4. Plot ---")
    plt.figure(figsize=(15, 6))
    
    plt.imshow(full_spectrogram, aspect='auto', origin='lower', cmap=colormap)
    
    plt.colorbar(label='Magnitude (dB)')
    plt.title(f"Full Signal Spectrogram (SNR = {target_snr} dB)")
    plt.xlabel("Time Step (approx)")
    plt.ylabel("Frequency Bin")
    
    if not os.path.exists('results'):
        os.makedirs('results')
        
    save_path = f'results/full_signal_stft_snr_{target_snr}dB.png'
    plt.tight_layout()
    plt.savefig(save_path)
    print(f"Plot saved to {save_path}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Visualize FHSS Signal STFT")
    parser.add_argument('--snr', type=int, default=10, help='SNR level in dB (default: 10)')
    parser.add_argument('--cmap', type=str, default='viridis', help='Matplotlib colormap (default: inferno)')
    
    args = parser.parse_args()
    visualize_full_signal(args.snr, args.cmap)
