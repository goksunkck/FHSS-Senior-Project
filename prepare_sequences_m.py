
import numpy as np
import h5py
import os
import sys

def main():
    print("--- Preparing M-SEQUENCE Sequences (Index-Based) ---")
    
    # Config
    dataset_path = os.path.join('data', 'synthetic', 'm_sequence_dataset.h5')
    output_path = os.path.join('data', 'synthetic', 'prepared_m_sequence_indices.h5')
    lookback_window = 12 
    
    if not os.path.exists(dataset_path):
        print(f"Error: {dataset_path} not found. Run generate_m_sequence_dataset.py first.")
        return

    with h5py.File(dataset_path, 'r') as f_in, h5py.File(output_path, 'w') as f_out:
        # Input
        X_all = f_in['X_data'] 
        Y_all = f_in['Y_data'] 
        SetIDs = f_in['Set_ID'] 
        
        # Meta
        try:
           num_hops_per_signal = int(f_in.attrs['numHopsPerSignal'])
        except:
           num_hops_per_signal = 256
           print("Warning: numHopsPerSignal missing, using 256")

        total_examples = X_all.shape[0]
        num_signals = total_examples // num_hops_per_signal
        
        print(f"Scanning {num_signals} signals for Train/Val/Test split...")
        
        # Arrays to hold start indices
        starts_train = []
        labels_train = []
        
        starts_val = []
        labels_val = []
        
        starts_test = []
        labels_test = []
        
        # Iterate Signal-by-Signal
        for sig_idx in range(num_signals):
            base = sig_idx * num_hops_per_signal
            
            # Check Set ID
            try:
                set_id = SetIDs[base][0]
            except:
                set_id = SetIDs[base]

            num_valid = num_hops_per_signal - lookback_window
            
            for i in range(num_valid):
                start_global = base + i
                target_global = base + i + lookback_window
                
                target_val = Y_all[target_global][0]
                
                if set_id == 0:
                    starts_train.append(start_global)
                    labels_train.append(target_val)
                elif set_id == 1:
                    starts_val.append(start_global)
                    labels_val.append(target_val)
                else:
                    starts_test.append(start_global)
                    labels_test.append(target_val)
            
            if (sig_idx+1) % 100 == 0:
                print(f"Scanned {sig_idx+1}/{num_signals}...")

        # Save to HDF5
        print(f"Saving indices... Train: {len(starts_train)}, Val: {len(starts_val)}, Test: {len(starts_test)}")
        
        f_out.create_dataset('sequence_starts_train', data=np.array(starts_train, dtype=int))
        f_out.create_dataset('YTrain', data=np.array(labels_train, dtype=int))
        
        f_out.create_dataset('sequence_starts_validation', data=np.array(starts_val, dtype=int))
        f_out.create_dataset('YValidation', data=np.array(labels_val, dtype=int))
        
        f_out.create_dataset('sequence_starts_test', data=np.array(starts_test, dtype=int))
        f_out.create_dataset('YTest', data=np.array(labels_test, dtype=int))
        
        # Also copy meta
        f_out.attrs['lookback_window'] = lookback_window
        
    print("Done.")

if __name__ == "__main__":
    main()
