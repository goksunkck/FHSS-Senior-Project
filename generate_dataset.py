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
    num_signals_per_snr = 100
    
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
        
        dset_X = f.create_dataset('X_data', (total_examples, 1024, 3), dtype='float32') # (N, Freq, Time)
        dset_Y = f.create_dataset('Y_data', (total_examples, 1), dtype='uint8')
        
        # Meta
        f.attrs['numHopsPerSignal'] = num_hops
        f.attrs['polynomial_degree'] = 10
        
        global_min = float('inf')
        global_max = float('-inf')
        
        write_idx = 0
        
        print(f"Generating {total_signals} signals...")
        
        for snr_db in snr_levels_db:
            print(f"  SNR: {snr_db} dB")
            
            for sig_idx in range(num_signals_per_snr):
                # Randomized Seed
                rand_seed = np.random.randint(0, 2, 10)
                while np.all(rand_seed == 0):
                    rand_seed = np.random.randint(0, 2, 10)
                
                # Gold Code
                gold_gen = GoldCodeGenerator(pn_poly1, pn_initial1, pn_poly2, rand_seed)
                # bits needed for hops: numHops * k
                hop_bits_needed = num_hops * k
                gold_seq = gold_gen.generate(hop_bits_needed) # (N, )
                
                # Map to Hops
                # reshape (k, numHops) or (numHops, k)?
                # MATLAB: reshape(seq, k, []) -> k rows. Transpose -> numHops rows.
                gold_seq_reshaped = gold_seq.reshape(k, -1, order='F').T # 'F' to match MATLAB reshape column-major
                # Actually MATLAB reshape fills columns first. 
                # goldCodeSequence is 1D vector.
                # reshape(v, k, []) -> 
                # col 1: v(1)..v(k)
                # col 2: v(k+1)..v(2k)
                # Transpose -> row 1: v(1)..v(k).
                # This corresponds to standard bits->int conversion.
                
                # In Python reshape: 'C' (default) fills rows first? No, last dim first.
                # 'F' fills first dim first (column-major).
                # So reshape(k, -1, order='F') is correct equivalent.
                
                # Convert bits to int
                # bi2de(..., 'left-msb') -> [1 0 0] = 4? OR [0 0 1] = 1?
                # MATLAB bi2de default is 'right-msb' (LSB first). 
                # Script used 'left-msb'. So first bit is MSB.
                # [b0 b1 b2] -> b0*4 + b1*2 + b2*1
                
                powers = 2 ** np.arange(k-1, -1, -1) # [4, 2, 1]
                hop_indices = gold_seq_reshaped.dot(powers).astype(int) # (numHops, )
                hop_freqs = hopset[hop_indices]
                
                # Data
                message_bits = np.random.randint(0, 2, num_bits)
                modulated_data = fsk_modulate(message_bits, sim_params)
                
                # Hops
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
                    
                    # STFT
                    stft_log = compute_stft_matrix(rx_signal, fs, stft_nfft, stft_win, stft_over)
                    
                    # Updates
                    g_min = np.min(stft_log)
                    g_max = np.max(stft_log)
                    if g_min < global_min: global_min = g_min
                    if g_max > global_max: global_max = g_max
                    
                    X_batch[h] = stft_log
                    Y_batch[h] = hop_indices[h]
                
                # Write batch
                dset_X[write_idx : write_idx+num_hops] = X_batch
                dset_Y[write_idx : write_idx+num_hops] = Y_batch
                write_idx += num_hops
                
        # Save Global Stats
        f.create_dataset('global_min', data=np.array([global_min]))
        f.create_dataset('global_max', data=np.array([global_max]))
        
        print(f"Done. Saved {write_idx} examples. Global Min: {global_min}, Global Max: {global_max}")

if __name__ == "__main__":
    main()
