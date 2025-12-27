import h5py
import numpy as np

# Check MATLAB dataset
matlab_file = 'data/synthetic/classification_dataset_stft_random_deg10.mat'

with h5py.File(matlab_file, 'r') as f:
    print("MATLAB Dataset Keys:", list(f.keys()))
    
    if 'X_data' in f:
        X = f['X_data']
        print(f"X_data shape: {X.shape}")
        print(f"X_data dtype: {X.dtype}")
        
        # Sample some data
        print("\nSampling 100 random frames...")
        indices = np.random.choice(X.shape[0], min(100, X.shape[0]), replace=False)
        indices.sort()
        
        # MATLAB stores as cell arrays, each cell contains a matrix
        # For HDF5, this is refs
        sample_values = []
        for idx in indices[:10]:  # Just first 10 to avoid too much loading
            try:
                ref = X[idx, 0]
                if isinstance(ref, h5py.h5r.Reference):
                    data = f[ref][:]
                    sample_values.append(data)
                    print(f"Frame {idx}: shape={data.shape}, min={data.min():.2f}, max={data.max():.2f}, mean={data.mean():.2f}")
            except Exception as e:
                print(f"Frame {idx}: Error - {e}")
        
        if sample_values:
            all_vals = np.concatenate([v.flatten() for v in sample_values])
            print(f"\nAggregate Stats (10 frames):")
            print(f"  Min: {all_vals.min():.2f}")
            print(f"  Max: {all_vals.max():.2f}")
            print(f"  Mean: {all_vals.mean():.2f}")
            print(f"  Std: {all_vals.std():.2f}")
    
    # Check for global stats
    if 'global_min' in f:
        print(f"\nGlobal Min (dataset): {f['global_min'][:]}")
    if 'global_max' in f:
        print(f"Global Max (dataset): {f['global_max'][:]}")
    
    # Check attributes
    print("\nAttributes:", dict(f.attrs))
