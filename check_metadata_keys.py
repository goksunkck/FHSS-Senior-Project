import h5py

# Check what keys are in the metadata file
metadata_file = 'data/synthetic/prepared_prediction_sequences_python.h5'

with h5py.File(metadata_file, 'r') as f:
    print("Keys in metadata file:")
    for key in f.keys():
        print(f"  {key}: shape={f[key].shape}, dtype={f[key].dtype}")
    
    print("\nAttributes:")
    for key in f.attrs.keys():
        print(f"  {key}: {f.attrs[key]}")
