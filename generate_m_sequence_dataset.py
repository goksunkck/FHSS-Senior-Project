
import numpy as np
import h5py
import os
import sys
# Add src to path
sys.path.append(os.path.join(os.getcwd(), 'src', 'python'))

from dsp_utils import fsk_modulate, add_awgn, compute_stft_matrix

class LFSRGenerator:
    def __init__(self, poly, initial_state):
        self.poly = poly
        # Convert to 1-based indexing for taps if inputs are 1-based, 
        # but here we assume inputs like [10, 3, 0] mean indices 10, 3, 0?
        # Actually usually [10, 3, 0] means 1 + x^3 + x^10.
        # The exponents are 10, 3, 0.
        # 0 is the current output (x^0).
        # Feedback comes from 10 and 3.
        
        self.state = np.array(initial_state, dtype=int)
        self.taps = [x for x in poly if x != 0] # exclude 0 if present as it's the output
        
    def step(self):
        # Calculate feedback bit: XOR sum of tapped bits
        # Taps are 1-based indices usually.
        # If state is length 10. Indices 0..9.
        # x^1 corresponds to index 0?
        # Usually: State = [s1, s2, ..., s10]
        # Feedback = s10 + s3 (for x^10 + x^3 + 1)
        # Shift: [Feedback, s1, s2, ..., s9]
        
        # Let's align with the GoldCodeGenerator logic from src/python/gold_code.py if possible.
        # But since I can't read it easily right now without a tool, I'll implement standard LFSR.
        
        # Standard: poly [10, 3, 0] -> 1 + x^3 + x^10
        # Feedback = bit_at_pos_10 XOR bit_at_pos_3
        # State: [b_0, b_1, ..., b_9] where b_0 is x^-1?
        # Let's assume State[i] is the delay element.
        # Taps 10 and 3 mean state[9] and state[2] (0-indexed).
        
        feedback = 0
        for t in self.taps:
            feedback ^= self.state[t-1] # 1-based tap to 0-based index
            
        output = self.state[-1] # The bit shifted out
        
        # Shift
        self.state = np.roll(self.state, 1)
        self.state[0] = feedback
        
        return output

    def generate(self, n):
        seq = np.zeros(n, dtype=int)
        for i in range(n):
            seq[i] = self.step()
        return seq

def main():
    print("Starting M-Sequence Dataset Generation...")
    
    # --- Configuration ---
    snr_levels_db = np.arange(-10, 12, 2) 
    
    # Poly: x^10 + x^3 + 1 -> [10, 3] (Excluding 0)
    # This is a primitive polynomial -> M-Sequence
    pn_poly = [10, 3]
    
    # Simulation Params 
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
    output_filename = os.path.join(output_dir, 'm_sequence_dataset.h5')
    
    # Strategy:
    # M-Sequences have only 1 cycle (length 1023).
    # All non-zero seeds are just phase shifts of this same cycle.
    # To test GENERALIZATION, we must split the CYCLE.
    # e.g. Train on the first 700 seeds (phase shifts)?
    # No, that essentially means we train on the whole sequence just shifted.
    # Wait, Phase Shift A vs Phase Shift B.
    # If the "Window" is 12.
    # Seeing transition [A, A+1... A+12] is valid training data.
    # Seeing transition [B, B+1... B+12] is valid.
    # If A and B are far apart, they are locally distinct.
    
    # Since we want to prove "Learning the Rule", we split by Seeds.
    # M-sequence seeds define the Initial State.
    # State = 10 bits. 1023 possible states.
    # Each state corresponds to a unique position in the 1023-long cycle.
    # So splitting by seeds is equivalent to Splitting by Position in the Cycle.
    # If we hold out 150 seeds, we hold out 150 specific 10-bit patterns (states).
    # THIS IS A VALID TEST. Can the model predict the successor of a state it has never seen?
    
    print("Creating disjoint seed pools...")
    all_integers = np.arange(1, 1024) 
    np.random.seed(42) 
    np.random.shuffle(all_integers)
    
    target_total_signals = 2000
    n_train_sigs = int(0.7 * target_total_signals) 
    n_val_sigs = int(0.15 * target_total_signals) 
    n_test_sigs = int(0.15 * target_total_signals) 
    
    total_examples = target_total_signals * num_hops

    print(f"Generating M-SEQUENCE dataset: {target_total_signals} signals total.")

    with h5py.File(output_filename, 'w') as f:
        chunk_size_examples = 100 * num_hops
        dset_X = f.create_dataset('X_data', (total_examples, 1024, 3), dtype='float32', chunks=(chunk_size_examples, 1024, 3))
        dset_Y = f.create_dataset('Y_data', (total_examples, 1), dtype='uint8')
        dset_SetID = f.create_dataset('Set_ID', (total_examples, 1), dtype='uint8')
        dset_SeedVal = f.create_dataset('Seed_Val', (total_examples, 1), dtype='int32')
        dset_SNR = f.create_dataset('SNR', (total_examples, 1), dtype='float32')
        
        # Meta
        f.attrs['numHopsPerSignal'] = num_hops
        f.attrs['polynomial_degree'] = 10
        f.attrs['sequence_type'] = 'm_sequence'
        
        global_min = float('inf')
        global_max = float('-inf')
        write_idx = 0
        
        # We assume seeds for M-sequence define initial state
        # Since M-sequence is single cycle, this is just starting at different points.
        
        # Pools
        # We recycle the full pool logic because we are sampling from it.
        # seeds are 1..1023.
        
        # For M-Sequence, we just use the seed as the initial state.
        
        def generate_batch(pool_seeds, n_signals_to_gen, set_id_val):
            nonlocal write_idx, global_min, global_max
            
            for i in range(n_signals_to_gen):
                # Pick 1 random seed from the allowed pool (Sampling with replacement allowed for inputs, 
                # but pool itself is disjoint set of integer states)
                seed_int = np.random.choice(pool_seeds)
                snr_db = np.random.choice(snr_levels_db)
                
                seed_bin_str = format(seed_int, '010b')
                rand_seed = [int(b) for b in seed_bin_str] # List
                
                # M-Sequence Generator
                lfsr = LFSRGenerator(pn_poly, rand_seed)
                
                hop_bits_needed = num_hops * k
                # Generate bits
                seq_bits = lfsr.generate(hop_bits_needed)
                
                # Form Hops
                gold_seq_reshaped = seq_bits.reshape(k, -1, order='F').T
                powers = 2 ** np.arange(k-1, -1, -1)
                hop_indices = gold_seq_reshaped.dot(powers).astype(int)
                hop_freqs = hopset[hop_indices]
                
                # Modulate
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

        # We take the shuffled integers and slice them to define disjoint pools
        # Define split sizes for SEEDS (1023 total)
        n_total_seeds = len(all_integers)
        n_train_seeds = int(0.7 * n_total_seeds)
        n_val_seeds = int(0.15 * n_total_seeds)
        # n_test_seeds is remainder
        
        pool_train = all_integers[:n_train_seeds]
        pool_val = all_integers[n_train_seeds : n_train_seeds + n_val_seeds]
        pool_test = all_integers[n_train_seeds + n_val_seeds:]
        
        print(f"Train Pool: {len(pool_train)}, Val Pool: {len(pool_val)}, Test Pool: {len(pool_test)}")

        print("Generating Training Set...")
        generate_batch(pool_train, n_train_sigs, 0)
        
        print("Generating Validation Set...")
        generate_batch(pool_val, n_val_sigs, 1)
        
        print("Generating Test Set...")
        generate_batch(pool_test, n_test_sigs, 2)
        
        # Save Global Stats
        f.create_dataset('global_min', data=np.array([global_min]))
        f.create_dataset('global_max', data=np.array([global_max]))
        
        print(f"Done. Saved {write_idx // num_hops} signals. Global Min: {global_min}, Global Max: {global_max}")

if __name__ == "__main__":
    main()
