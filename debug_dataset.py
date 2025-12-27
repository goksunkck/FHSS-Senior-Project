
import h5py
import numpy as np

def check_consistency():
    prep_file = 'data/synthetic/prepared_prediction_sequences.mat'
    data_file = 'data/synthetic/classification_dataset_stft_awgn.mat'
    
    print(f"--- Checking {prep_file} (HDF5 mode) ---")
    try:
        with h5py.File(prep_file, 'r') as f:
            print("Keys:", list(f.keys()))
            if 'sequence_starts_train' in f:
                # MATLAB saves vectors as (N, 1) or (1, N).
                seq_starts = f['sequence_starts_train'][()]
                print(f"sequence_starts_train shape: {seq_starts.shape}")
                
                # Check max index
                # Note: MATLAB indices are 1-based double.
                max_idx = np.max(seq_starts)
                min_idx = np.min(seq_starts)
                print(f"Validation of Indices (1-based from MATLAB):")
                print(f"  Min Start Index: {min_idx}")
                print(f"  Max Start Index: {max_idx}")
                
                # Store for comparison
                max_required_idx = max_idx # This is the row index needed in X_data
            else:
                print("Could not find sequence_starts_train")
                return
    except Exception as e:
        print(f"Error reading prep file: {e}")
        return

    print(f"\n--- Checking {data_file} (HDF5 mode) ---")
    try:
        with h5py.File(data_file, 'r') as f:
            print("Keys:", list(f.keys()))
            if 'X_data' in f:
                x_dset = f['X_data']
                print(f"X_data shape: {x_dset.shape}")
                print(f"X_data type: {x_dset.dtype}")
                
                total_rows = x_dset.shape[0] if len(x_dset.shape) > 0 else 0
                
                print(f"\n--- Consistency Check ---")
                print(f"Prepared Metadata asks for index up to: {max_required_idx}")
                print(f"Actual Dataset has rows: {total_rows}")
                
                if max_required_idx > total_rows:
                    print("!!! CRITICAL MISMATCH !!!")
                    print("The metadata points to indices that define logical positions larger than the dataset.")
                    print("Conclusion: 'prepared_prediction_sequences.mat' is OUT OF SYNC with 'classification_dataset_stft_awgn.mat'.")
                    print("Action: You must run 'prepare_prediction_sequences.m' in MATLAB to regenerate the metadata.")
                else:
                    print("Indices are logically valid.")
            else:
                print("X_data key not found!")
    except Exception as e:
        print(f"Error reading data file: {e}")

if __name__ == "__main__":
    check_consistency()
