import numpy as np
import h5py
import os

def main():
    print("--- Preparing Sequences (Python) ---")
    
    # Config
    dataset_path = os.path.join('data', 'synthetic', 'classification_dataset_stft_random_deg10_python.h5')
    output_path = os.path.join('data', 'synthetic', 'prepared_prediction_sequences_python.h5')
    lookback_window = 35
    
    if not os.path.exists(dataset_path):
        print(f"Error: {dataset_path} not found.")
        return

    with h5py.File(dataset_path, 'r') as f_in:
        # Read Meta
        num_hops_per_signal = int(f_in.attrs['numHopsPerSignal'])
        # Read Y data for splitting
        Y_data = f_in['Y_data'][:] # Load all labels to memory (small enough: 330 * 256 bytes)
        total_examples = Y_data.shape[0]
        
        num_signals = total_examples // num_hops_per_signal
        
        print(f"Loaded {total_examples} examples ({num_signals} signals).")
        
        # Sequences
        num_seq_per_sig = num_hops_per_signal - lookback_window
        total_sequences = num_signals * num_seq_per_sig
        
        sequence_starts = np.zeros(total_sequences, dtype=np.int32)
        Y_labels = np.zeros(total_sequences, dtype=np.uint8)
        signal_ids = np.zeros(total_sequences, dtype=np.int32)
        
        seq_idx = 0
        for i in range(num_signals):
            start_idx = i * num_hops_per_signal
            try:
                sig_y = Y_data[start_idx : start_idx + num_hops_per_signal].flatten()
            except:
                sig_y = Y_data[start_idx : start_idx + num_hops_per_signal, 0]
                
            for j in range(num_seq_per_sig):
                # sequence start index (in X_data)
                # MATLAB was 1-based, Python 0-based.
                # X_data[seq_start : seq_start + lookback] is the input
                # Target is Y[seq_start + lookback]
                
                # Input: [t, t+1, ..., t+LB-1] -> Predict Y[t+LB] ?
                # prepare_randomized_sequences.m: 
                # sequence_starts(seq_counter) = signal_start_idx + j - 1; (j=1..N)
                # Y_labels(seq_counter) = current_Y_sequence(j + lookback_window);
                
                # Equivalent in Python:
                # Sequence j (0..N-1):
                # Start index in X: signal_start + j
                # Label index: signal_start + j + lookback
                
                sequence_starts[seq_idx] = start_idx + j
                Y_labels[seq_idx] = sig_y[j + lookback_window]
                signal_ids[seq_idx] = i
                seq_idx += 1
                
        # Split (Signal-wise)
        np.random.seed(42)
        perm = np.random.permutation(num_signals)
        num_train = int(0.7 * num_signals)
        num_val = int(0.15 * num_signals)
        
        train_sigs = perm[:num_train]
        val_sigs = perm[num_train : num_train + num_val]
        test_sigs = perm[num_train + num_val:]
        
        # Create masks
        split_map = np.zeros(num_signals, dtype=np.uint8) # 0:Train, 1:Val, 2:Test
        split_map[val_sigs] = 1
        split_map[test_sigs] = 2
        
        seq_splits = split_map[signal_ids]
        
        idx_train = seq_splits == 0
        idx_val = seq_splits == 1
        idx_test = seq_splits == 2
        
        print(f"Split: {np.sum(idx_train)} Train, {np.sum(idx_val)} Val, {np.sum(idx_test)} Test sequences.")
        
        # Save
        with h5py.File(output_path, 'w') as f_out:
            f_out.create_dataset('sequence_starts_train', data=sequence_starts[idx_train])
            f_out.create_dataset('YTrain', data=Y_labels[idx_train])
            
            f_out.create_dataset('sequence_starts_validation', data=sequence_starts[idx_val])
            f_out.create_dataset('YValidation', data=Y_labels[idx_val])
            
            f_out.create_dataset('sequence_starts_test', data=sequence_starts[idx_test])
            f_out.create_dataset('YTest', data=Y_labels[idx_test])
            
            # Copy stats
            if 'global_min' in f_in:
                f_out.create_dataset('global_min', data=f_in['global_min'][:])
                f_out.create_dataset('global_max', data=f_in['global_max'][:])
            
            # Save meta
            f_out.attrs['dataset_filename'] = dataset_path
            f_out.attrs['lookback_window'] = lookback_window
            
            # Class values
            classes = np.unique(Y_data)
            f_out.create_dataset('class_values', data=classes)
            
    print(f"Saved prepared sequences to {output_path}")

if __name__ == "__main__":
    main()
