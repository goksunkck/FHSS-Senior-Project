import h5py
import numpy as np
import matplotlib.pyplot as plt

def verify_data():
    dataset_path = 'data/synthetic/prepared_prediction_sequences_python.h5'
    source_path = 'data/synthetic/classification_dataset_stft_random_deg10_python.h5'
    
    print(f"--- Verifying Data Integrity ---")
    print(f"Dataset: {dataset_path}")
    
    with h5py.File(dataset_path, 'r') as f, h5py.File(source_path, 'r') as f_src:
        # Load sample indices
        starts = f['sequence_starts_train'][:100]
        labels = f['YTrain'][:100]
        
        # Load full source data to check alignment
        X_full = f_src['X_data']
        Y_full = f_src['Y_data']
        
        lookback = f.attrs['lookback_window']
        print(f"Lookback Window: {lookback}")
        
        correct_alignments = 0
        total_checks = 10
        
        print("\nChecking first 10 training sequences...")
        for i in range(total_checks):
            start_idx = starts[i]
            target_y = labels[i]
            
            # Reconstruct what the input sequence labels WERE
            # We don't have X labels directly in 'prepared' file, but we can look up Y_full
            # The input sequence corresponds to Y_full[start_idx : start_idx + lookback]
            # The target should be Y_full[start_idx + lookback]
            
            input_y_seq = Y_full[start_idx : start_idx + lookback].flatten()
            actual_next_step = Y_full[start_idx + lookback][0]
            
            print(f"Seq {i}: Input classes: {input_y_seq[-5:]} ... -> Target: {target_y}")
            print(f"       Actual Next in Src: {actual_next_step}")
            
            if target_y == actual_next_step:
                correct_alignments += 1
            else:
                print(f"   !!! MISMATCH !!! Target {target_y} != Actual {actual_next_step}")
        
        print(f"\nAlignment Check: {correct_alignments}/{total_checks} correct.")
        
        # Check Value Ranges
        print("\nStats:")
        print(f"YTrain Unique: {np.unique(f['YTrain'][:])}")
        
if __name__ == "__main__":
    verify_data()
