import torch
import matplotlib.pyplot as plt
import numpy as np
from train_frequency_predictor_colab import RAMDataset, load_mat_file
import h5py
import os

def debug_loader():
    print("--- Debugging DataLoader ---")
    
    # 1. Load Metadata (Copied from train script)
    prep_file_python = 'data/synthetic/prepared_prediction_sequences_python.h5'
    
    if not os.path.exists(prep_file_python):
        print("Error: Metadata file not found.")
        return

    prep_data = {}
    with h5py.File(prep_file_python, 'r') as f:
        prep_data['sequence_starts_train'] = f['sequence_starts_train'][:]
        prep_data['YTrain'] = f['YTrain'][:]
        # Attributes
        prep_data['lookback_window'] = int(f.attrs['lookback_window']) if 'lookback_window' in f.attrs else int(f['lookback_window'][0])
        
        dset_name = f.attrs['dataset_filename'] if 'dataset_filename' in f.attrs else f['dataset_filename'][:]
        if isinstance(dset_name, (bytes, np.bytes_)):
             prep_data['dataset_filename'] = dset_name.decode('utf-8')
        elif not isinstance(dset_name, str):
             prep_data['dataset_filename'] = "".join([chr(c) for c in dset_name.flatten()])
        else:
             prep_data['dataset_filename'] = str(dset_name)
             
        prep_data['class_values'] = f['class_values'][:]

    # 2. Get Norm Stats
    dataset_filename = prep_data['dataset_filename']
    print(f"Dataset: {dataset_filename}")
    
    g_min = -50
    g_max = 50
    with h5py.File(dataset_filename, 'r') as f:
        if 'global_min' in f.attrs:
             g_min = float(f.attrs['global_min'])
             g_max = float(f.attrs['global_max'])
        elif 'global_min' in f:
             g_min = float(f['global_min'][0])
             g_max = float(f['global_max'][0])
    
    print(f"Norm: Min={g_min}, Max={g_max}")
    g_range = g_max - g_min
    
    # 3. Create Dataset
    # Map labels just for train
    class_values = prep_data['class_values']
    class_map = {val: i for i, val in enumerate(class_values)}
    YTrain_mapped = np.array([class_map[y] for y in prep_data['YTrain']], dtype=np.int64)
    
    # Subset for speed (first 100)
    dataset = RAMDataset(
        dataset_filename, 
        prep_data['sequence_starts_train'][:100], 
        YTrain_mapped[:100], 
        prep_data['lookback_window'], 
        g_min, 
        g_range
    )
    
    # 4. Get one item and visualize
    print("Fetching index 0...")
    seq, label = dataset[0]
    # Seq shape: (35, 1, 1024, 3)
    
    print(f"Sequence Shape: {seq.shape}")
    print(f"Label: {label} (Original: {prep_data['YTrain'][0]})")
    print(f"Data Range: Min={seq.min()}, Max={seq.max()}, Mean={seq.mean()}")
    
    # Check if data is all zeros or NaNs
    if torch.isnan(seq).any():
        print("WARNING: NaNs found in sequence!")
    if seq.max() == seq.min():
        print("WARNING: Sequence is constant (flat)!")
        
    # Visualize frame 0
    # (1, 1024, 3) -> (1024, 3)
    frame0 = seq[0, 0, :, :].numpy()
    
    plt.figure(figsize=(10, 5))
    plt.subplot(1, 2, 1)
    plt.imshow(frame0, aspect='auto', origin='lower')
    plt.title(f"Frame 0 (Normalized)")
    plt.colorbar()
    
    # Visualize frame 10
    frame10 = seq[10, 0, :, :].numpy()
    plt.subplot(1, 2, 2)
    plt.imshow(frame10, aspect='auto', origin='lower')
    plt.title(f"Frame 10 (Normalized)")
    plt.colorbar()
    
    plt.savefig('results/debug_dataloader_sample.png')
    print("Saved results/debug_dataloader_sample.png")

if __name__ == "__main__":
    debug_loader()
