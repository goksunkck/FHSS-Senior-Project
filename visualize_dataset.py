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
    
    # Check for Python version first
    prep_file_python = 'data/synthetic/prepared_prediction_sequences_python.h5'
    prep_file_mat = 'data/synthetic/prepared_prediction_sequences.mat'
    
    prep_file = prep_file_python if os.path.exists(prep_file_python) else prep_file_mat
    print(f"Loading metadata from: {prep_file}")

    dataset_filename = None
    sequence_starts_train = None
    YTrain = None
    
    if prep_file.endswith('.h5'):
        # Load HDF5 directly
        with h5py.File(prep_file, 'r') as f:
            sequence_starts_train = f['sequence_starts_train'][:]
            YTrain = f['YTrain'][:]
            dataset_filename = f.attrs['dataset_filename']
            # Decode bytes if needed
            if isinstance(dataset_filename, bytes):
                dataset_filename = dataset_filename.decode('utf-8')
    else:
        # Legacy MAT load
        try:
            prep_data = load_mat_file(prep_file, 'prep_data')
            sequence_starts_train = prep_data['sequence_starts_train'].flatten()
            YTrain = prep_data['YTrain'].flatten()
            dataset_filename = prep_data['dataset_filename']
            if isinstance(dataset_filename, np.ndarray):
                dataset_filename = "".join(map(chr, dataset_filename.flatten()))
        except Exception as e:
            print(f"Error loading legacy MAT file: {e}")
            return
            
    print(f"Dataset filename referenced: {dataset_filename}")

    print("--- 2. Load base STFT dataset into memory ---")
    source_dataset = dataset_filename
    
    if not os.path.exists(source_dataset):
        print(f"Error: Source dataset not found at {source_dataset}.")
        return

    # Check file type of source
    is_python_source = source_dataset.endswith('.h5') and 'python' in source_dataset
    
    hf_source = None
    X_data_source = None
    X_data_in_memory = None
    is_direct = False

    if is_python_source:
         print("Source is Python HDF5. Opening file handle...")
         hf_source = h5py.File(source_dataset, 'r')
         X_data_source = hf_source['X_data']
         is_direct = True
    else:
        try:
            data_contents = load_mat_file(source_dataset, 'X_data')
            X_data_in_memory = data_contents['X_data']
            print("Dataset loaded into memory.")
            is_direct = False

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
                
        except Exception as e:
            print(f"Error loading source dataset: {e}")
            return

    # --- 3. Visualize Random Samples ---
    print("--- 3. Visualizing Random Samples ---")
    
    try:
        num_samples = 9
        indices = random.sample(range(len(sequence_starts_train)), num_samples)
        
        plt.figure(figsize=(15, 15))
        
        for i, idx in enumerate(indices):
            # sequence_starts_train contains the START index of a sequence in X_data
            start_idx = int(sequence_starts_train[idx])
            
            # Adjust for MATLAB 1-based indexing if NOT using Python format
            if not is_direct:
                start_idx = start_idx - 1
                if start_idx < 0: start_idx = 0
                
                # Get the image
                img_cell = X_data_in_memory[0, start_idx] # Assuming (1, N)
                
                if isinstance(img_cell, np.ndarray):
                    img = img_cell
                else:
                    img = img_cell[0]
            else:
                # Direct Access from HDF5
                img = X_data_source[start_idx]
                
            plot_img = np.array(img)
            
            # If 3D (e.g. 1024, 3, 1), squeeze
            if plot_img.ndim == 3:
                 plot_img = np.squeeze(plot_img)
                 
            # Now we likely have (1024, 3).
            # Check if transposed
            if plot_img.shape[0] < plot_img.shape[1]:
                 plot_img = plot_img.T
                 
            # Normalize for display if needed (0-1)
            p_min, p_max = plot_img.min(), plot_img.max()
            if p_max > p_min:
                plot_img = (plot_img - p_min) / (p_max - p_min)
                
            plt.subplot(3, 3, i + 1)
            plt.imshow(plot_img, aspect='auto', cmap='viridis')
            
            label_val = YTrain[idx]
            plt.title(f"Label: {label_val}\nShape: {img.shape}")
            plt.axis('off')

        if not os.path.exists('results'):
            os.makedirs('results')
            
        save_path = 'results/dataset_sample.png'
        plt.tight_layout()
        plt.savefig(save_path)
        print(f"Visualization saved to {save_path}")

    finally:
        if hf_source:
            hf_source.close()

if __name__ == "__main__":
    visualize()
