import numpy as np
import h5py
import os
import sys
# Add src to path
sys.path.append(os.path.join(os.getcwd(), 'src', 'python'))

from gold_code import GoldCodeGenerator
from dsp_utils import fsk_modulate, add_awgn, compute_stft_matrix

def main():
    print("Starting Python Dataset Generation...")
    
    # --- Configuration ---
    snr_levels_db = np.arange(-10, 12, 2) # -10 to 10 step 2
    num_signals_per_snr = 100 ### NOT USED
    
    # Poly1: x^10 + x^3 + 1 -> [10, 3, 0]
    pn_poly1 = [10, 3, 0]
    pn_initial1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 1]
    
    # Poly2: x^10 + x^7 + 1 -> [10, 7, 0]
    pn_poly2 = [10, 7, 0]
    
    # Simulation Params (matching MATLAB)
    fs = 10e6
    num_bits = 512
    M = 2
    symbol_rate = fs * 0.01
    freq_separation = symbol_rate
    samples_per_symbol = int(round(fs / symbol_rate))
    k = 3
    num_hops = 256
    num_channels = 2**k
    spacing = 2 * freq_separation
    base_freq = 2e6
    hopset = np.arange(num_channels) * spacing + base_freq
    
    sim_params = {
        'fs': fs,
        'M': M,
        'freqSeparation': freq_separation,
        'samplesPerSymbol': samples_per_symbol
    }
    
    # STFT Params
    stft_win = 128
    stft_over = 96
    stft_nfft = 1024
    
    samples_per_hop = int(np.floor((num_bits / np.log2(M) * samples_per_symbol) / num_hops))
    
    # Output file
    output_dir = os.path.join('data', 'synthetic')
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    output_filename = os.path.join(output_dir, 'classification_dataset_stft_random_deg10_python.h5')
    
    total_signals = len(snr_levels_db) * num_signals_per_snr
    total_examples = total_signals * num_hops
    
    # HDF5 Creation
    # We'll save as (Total, 1024, 3) directly.
    # STFT is complex? MATLAB script saved LogMag (real).
    # We stick to LogMag.
    
    with h5py.File(output_filename, 'w') as f:
        # Create datasets
        # Using a flattened shape for X like (Total, 1024, 3)
        # Assuming STFT output is 1024 x 3
        
        # NOTE: We need to confirm STFT time dimension. 
        # sample_per_hop size?
        # sample_per_hop = floor(512 / 1 * 100 / 256) = floor(200) = 200 samples?
        # fs=10e6, symRate=1e5 -> 100 samples/sym.
        # 512 bits -> 512 symbols -> 51200 samples.
        # 256 hops -> 200 samples/hop.
        # Window=128, Overlap=96.
        # L = 200.
        # 1st window: 1..128.
        # 2nd window: (128-96)+1 = 33..160.
        # 3rd window: (160-96)+1 = 65..192.
        # 4th window: 97..224 (Out of bounds).
        # So typically 3 columns. Correct.
        
    # --- POOL BASED GENERATION (Strictly Disjoint) ---
        # 1. Generate all possible non-zero seeds (1 to 1023)
        # 2. Shuffle them
        # 3. Partition into Train, Val, Test sets
        # 4. Generate data based on these sets
        
        print("Creating disjoint seed pools...")
        all_integers = np.arange(1, 1024) # 1 to 1023 (1023 seeds)
        np.random.seed(42) # Fixed seed for reproducibility of splits
        np.random.shuffle(all_integers)
        
        # Split: 70% Train, 15% Val, 15% Test
        n_total = len(all_integers)
        n_train = int(0.7 * n_total) # ~716
        n_val = int(0.15 * n_total)  # ~153
        n_test = n_total - n_train - n_val # ~154
        
        seeds_train = all_integers[:n_train]
        seeds_val = all_integers[n_train:n_train+n_val]
        seeds_test = all_integers[n_train+n_val:]
        
        print(f"Pool sizes: Train={len(seeds_train)}, Val={len(seeds_val)}, Test={len(seeds_test)}")
        
        # Strategy: Random Sampling from Disjoint Pools
        # Target: ~2000 signals total (manageable RAM size ~6GB when loaded)
        
        target_total_signals = 2000
        n_train_sigs = int(0.7 * target_total_signals) # 1400
        n_val_sigs = int(0.15 * target_total_signals)  # 300
        n_test_sigs = int(0.15 * target_total_signals) # 300
        
        # Adjust total to sum
        target_total_signals = n_train_sigs + n_val_sigs + n_test_sigs
        total_examples = target_total_signals * num_hops

        print(f"Generating OPTIMIZED dataset: {target_total_signals} signals total.")
        print(f"  Train: {n_train_sigs} | Val: {n_val_sigs} | Test: {n_test_sigs}")

        # Dataset Creation
        chunk_size_examples = 100 * num_hops
        dset_X = f.create_dataset('X_data', (total_examples, 1024, 3), dtype='float32', chunks=(chunk_size_examples, 1024, 3))
        dset_Y = f.create_dataset('Y_data', (total_examples, 1), dtype='uint8')
        dset_SetID = f.create_dataset('Set_ID', (total_examples, 1), dtype='uint8')
        dset_SeedVal = f.create_dataset('Seed_Val', (total_examples, 1), dtype='int32')
        dset_SNR = f.create_dataset('SNR', (total_examples, 1), dtype='float32')
        
        # Meta - IMPORTANT: Re-adding these so prepare_sequences can read them
        f.attrs['numHopsPerSignal'] = num_hops
        f.attrs['polynomial_degree'] = 10
        
        global_min = float('inf')
        global_max = float('-inf')
        write_idx = 0

        def generate_batch(pool_seeds, n_signals_to_gen, set_id_val):
            nonlocal write_idx, global_min, global_max
            
            for i in range(n_signals_to_gen):
                # Pick 1 random seed from the allowed pool
                seed_int = np.random.choice(pool_seeds)
                # Pick 1 random SNR
                snr_db = np.random.choice(snr_levels_db)
                
                # Generate Gold Code
                seed_bin_str = format(seed_int, '010b')
                rand_seed = np.array([int(b) for b in seed_bin_str])
                gold_gen = GoldCodeGenerator(pn_poly1, pn_initial1, pn_poly2, rand_seed)
                
                hop_bits_needed = num_hops * k
                gold_seq = gold_gen.generate(hop_bits_needed)
                
                # Map to Hops
                gold_seq_reshaped = gold_seq.reshape(k, -1, order='F').T
                powers = 2 ** np.arange(k-1, -1, -1)
                hop_indices = gold_seq_reshaped.dot(powers).astype(int)
                hop_freqs = hopset[hop_indices]
                
                # Modulate & Channel
                message_bits = np.random.randint(0, 2, num_bits)
                modulated_data = fsk_modulate(message_bits, sim_params)
                
                X_batch = np.zeros((num_hops, 1024, 3), dtype=np.float32)
                Y_batch = np.zeros((num_hops, 1), dtype=np.uint8)
                t_vec = np.arange(samples_per_hop) / fs
                
                for h in range(num_hops):
                    start = h * samples_per_hop
                    end = (h + 1) * samples_per_hop
                    segment = modulated_data[start:end]
                    carrier = np.exp(1j * 2 * np.pi * hop_freqs[h] * t_vec)
                    hop_signal = segment * carrier
                    rx_signal = add_awgn(hop_signal, snr_db)
                    
                    stft_log = compute_stft_matrix(rx_signal, fs, stft_nfft, stft_win, stft_over)
                    
                    # Update Min/Max
                    g_min = np.min(stft_log)
                    g_max = np.max(stft_log)
                    if g_min < global_min: global_min = g_min
                    if g_max > global_max: global_max = g_max
                    
                    X_batch[h] = stft_log
                    Y_batch[h] = hop_indices[h]

                # Write
                end_idx = write_idx + num_hops
                dset_X[write_idx : end_idx] = X_batch
                dset_Y[write_idx : end_idx] = Y_batch
                dset_SetID[write_idx : end_idx] = set_id_val
                dset_SeedVal[write_idx : end_idx] = seed_int 
                dset_SNR[write_idx : end_idx] = snr_db
                write_idx = end_idx
                
                if (i+1) % 50 == 0:
                     print(f"  [Set {set_id_val}] Generated {i+1}/{n_signals_to_gen} signals...")

        print("Generating Training Set...")
        generate_batch(seeds_train, n_train_sigs, 0)
        
        print("Generating Validation Set...")
        generate_batch(seeds_val, n_val_sigs, 1)
        
        print("Generating Test Set...")
        generate_batch(seeds_test, n_test_sigs, 2)
        
        # Save Global Stats
        f.create_dataset('global_min', data=np.array([global_min]))
        f.create_dataset('global_max', data=np.array([global_max]))
        
        print(f"Done. Saved {write_idx // num_hops} signals. Global Min: {global_min}, Global Max: {global_max}")

if __name__ == "__main__":
    main()
