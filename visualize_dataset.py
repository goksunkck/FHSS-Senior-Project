import numpy as np
import scipy.io
import h5py
import matplotlib.pyplot as plt
import os
import random

def load_mat_file(filepath, variable_name):
    """Loads a specific variable from a .mat file."""
    try:
        return scipy.io.loadmat(filepath)
    except NotImplementedError:
        # If scipy.io.loadmat fails, it might be a v7.3 .mat file.
        # h5py is needed for this.
        with h5py.File(filepath, 'r') as f:
            data = {}
            for k, v in f.items():
                if isinstance(v, h5py.Dataset): # Check if it's a dataset
                    # Dereference and copy data to make it mutable
                    value = v[()]
                    # Handle MATLAB's tendency to save everything as arrays
                    if isinstance(value, np.ndarray):
                        if value.shape == (1, 1):
                            data[k] = value[0, 0] # Unpack scalar
                        else:
                            data[k] = value
                    else:
                        data[k] = value
            return data

def visualize():
    print("--- 1. Load prepared data ---")
    prep_file = 'data/synthetic/prepared_prediction_sequences.mat'
    try:
        prep_data = load_mat_file(prep_file, 'prep_data')
    except FileNotFoundError:
        print(f"Error: Prepared dataset not found at {prep_file}.")
        return

    # Extract needed variables
    # Note: Keys might vary slightly depending on how they were saved/loaded
    # We try to be robust
    
    # Helper to safely get data whether it's in a struct or direct
    def get_val(data, key):
        if key in data:
            return data[key]
        # If prep_data was loaded as a dict containing 'prep_data' struct (scipy behavior sometimes)
        if 'prep_data' in data:
             # This part is tricky without knowing exact structure, but let's assume flat load for now
             # based on previous script
             pass
        return None

    # Based on train_frequency_predictor.py, prep_data seems to be the dictionary itself if loaded via h5py logic
    # or scipy.io.loadmat returns a dict.
    
    # Let's assume the keys are directly available or inside 'prep_data' key if scipy loaded it.
    # But the previous script used `prep_data['sequence_starts_train']`, implying `prep_data` IS the dictionary.
    
    sequence_starts_train = prep_data.get('sequence_starts_train')
    if sequence_starts_train is None:
        # Try looking inside a 'prep_data' key if it exists (scipy.io.loadmat behavior)
        if 'prep_data' in prep_data:
             # This would be a structured array
             # It's safer to rely on the previous script's working assumption:
             # prep_data IS the container.
             pass

    # Re-using logic from train_frequency_predictor.py which seemed to work
    sequence_starts_train = prep_data['sequence_starts_train'].flatten()
    YTrain = prep_data['YTrain'].flatten()
    
    dataset_filename = prep_data['dataset_filename']
    if isinstance(dataset_filename, np.ndarray):
        dataset_filename = "".join(map(chr, dataset_filename.flatten()))
    
    class_values = prep_data['class_values'].flatten()
    
    print(f"Dataset filename referenced: {dataset_filename}")

    print("--- 2. Load base STFT dataset into memory ---")
    source_dataset = dataset_filename
    try:
        data_contents = load_mat_file(source_dataset, 'X_data')
        X_data_in_memory = data_contents['X_data']
        print("Dataset loaded into memory.")

        # Resolve HDF5 references if needed
        if hasattr(X_data_in_memory, 'dtype') and X_data_in_memory.dtype == h5py.ref_dtype:
            print("Resolving HDF5 object references for X_data...")
            resolved_x_data = np.empty(X_data_in_memory.shape, dtype=object)
            with h5py.File(source_dataset, 'r') as f:
                for i in range(X_data_in_memory.size):
                    ref = X_data_in_memory.flat[i]
                    resolved_x_data.flat[i] = f[ref][()]
            X_data_in_memory = resolved_x_data
            print("References resolved.")
            
    except FileNotFoundError:
        print(f"Error: Source dataset not found at {source_dataset}.")
        return

    # --- 3. Visualize Random Samples ---
    print("--- 3. Visualizing Random Samples ---")
    
    num_samples = 9
    indices = random.sample(range(len(sequence_starts_train)), num_samples)
    
    plt.figure(figsize=(15, 15))
    
    for i, idx in enumerate(indices):
        # sequence_starts_train contains the START index of a sequence in X_data
        # Let's visualize the FIRST frame of that sequence
        start_idx_matlab = sequence_starts_train[idx]
        start_idx_python = start_idx_matlab - 1
        
        # Get the image
        img_cell = X_data_in_memory[0, start_idx_python]
        
        if isinstance(img_cell, np.ndarray):
            img = img_cell
        else:
            img = img_cell[0]
            
        # img is likely (3, 1024) or (3, 1024, 1) or similar.
        # We want to visualize it.
        # If it's (3, 1024), it's 3 channels (RGB?) x 1024 time/freq bins?
        # Actually, usually STFT is (Freq, Time).
        # Let's check shape.
        
        # From train script: "Transpose from (3, 1024, 1) to (1024, 3, 1)"
        # So original is likely (3, 1024) or (3, 1024, 1).
        # If we want to plot it as an image, we might want (1024, 3) or just plot one channel.
        
        if img.ndim == 3:
             # (H, W, C) or (C, H, W)?
             # Train script: "Transpose from (3, 1024, 1) to (1024, 3, 1)"
             # So it seems original is (3, 1024, 1) -> (Channels, Width?, 1)
             # Let's assume (Channels, Time)
             pass
        
        # Let's just try to plot it.
        # If it has 3 channels, we can try plotting as RGB if dimensions match.
        # If it's (3, 1024), it's very thin.
        
        # Let's print shape for debugging in the title
        shape_str = str(img.shape)
        
        # Prepare for plotting
        # If (3, 1024), let's transpose to (1024, 3) for better vertical view?
        # Or maybe it's (Freq, Time) = (1024, 3)? No, 3 is too small for time.
        # Maybe 3 is channels (Real, Imag, Mag)?
        
        plot_img = img
        if img.ndim == 3 and img.shape[2] == 1:
            plot_img = img[:, :, 0] # Drop singleton
            
        # If (3, 1024), transpose to (1024, 3)
        if plot_img.shape[0] == 3 and plot_img.shape[1] > 3:
            plot_img = plot_img.T
            
        # Normalize for display if needed (0-1)
        p_min, p_max = plot_img.min(), plot_img.max()
        if p_max > p_min:
            plot_img = (plot_img - p_min) / (p_max - p_min)
            
        plt.subplot(3, 3, i + 1)
        plt.imshow(plot_img, aspect='auto', cmap='viridis')
        
        label_val = YTrain[idx]
        plt.title(f"Label: {label_val}\nShape: {shape_str}")
        plt.axis('off')

    if not os.path.exists('results'):
        os.makedirs('results')
        
    save_path = 'results/dataset_sample.png'
    plt.tight_layout()
    plt.savefig(save_path)
    print(f"Visualization saved to {save_path}")

if __name__ == "__main__":
    visualize()
