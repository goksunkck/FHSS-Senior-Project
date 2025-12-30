"""
Test Berlekamp-Massey Algorithm on Real FHSS Dataset

This script:
1. Loads a signal from your dataset
2. Observes the first N frequency hops
3. Reconstructs the Gold code LFSR using Berlekamp-Massey
4. Predicts future hops
5. Compares predictions with actual hops
"""

import h5py
import numpy as np
import sys
import os

sys.path.append(os.path.join(os.getcwd(), 'src', 'python'))

from berlekamp_massey import reconstruct_gold_code_lfsr
from gold_code import LFSR, GoldCodeGenerator

def test_on_dataset():
    dataset_file = 'data/synthetic/classification_dataset_stft_random_deg10_python.h5'
    
    with h5py.File(dataset_file, 'r') as f:
        Y_data = f['Y_data']
        
        # Get one complete signal (256 hops)
        signal_hops = Y_data[0:256, 0]
        
        print("="*70)
        print("Testing Berlekamp-Massey on Real FHSS Dataset")
        print("="*70)
        
        # Test with different observation lengths
        for num_observed in [20, 35, 50]:
            print(f"\n{'='*70}")
            print(f"Observing first {num_observed} hops")
            print(f"{'='*70}")
            
            observed_hops = signal_hops[:num_observed]
            actual_future_hops = signal_hops[num_observed:num_observed+20]
            
            print(f"\nObserved hop sequence: {observed_hops}")
            
            # Reconstruct LFSR
            try:
                polynomial, initial_state, degree = reconstruct_gold_code_lfsr(observed_hops, k=3)
                
                # Now predict future hops
                print(f"\nPredicting next 20 hops using reconstructed LFSR...")
                
                # Generate bit sequence
                lfsr = LFSR(polynomial, initial_state)
                
                # Skip the observed bits
                num_observed_bits = num_observed * 3
                _ = lfsr.generate(num_observed_bits)
                
                # Generate next bits
                future_bits = lfsr.generate(20 * 3)  # 20 hops × 3 bits/hop
                
                # Convert to hop indices
                predicted_hops = []
                for i in range(0, len(future_bits), 3):
                    hop_bits = future_bits[i:i+3]
                    hop_index = hop_bits[0] * 4 + hop_bits[1] * 2 + hop_bits[2]
                    predicted_hops.append(hop_index)
                
                predicted_hops = np.array(predicted_hops)
                
                # Compare
                print(f"\nActual next 20 hops:    {actual_future_hops}")
                print(f"Predicted next 20 hops: {predicted_hops}")
                
                matches = np.sum(predicted_hops == actual_future_hops)
                accuracy = matches / len(actual_future_hops) * 100
                
                print(f"\nPrediction Accuracy: {matches}/{len(actual_future_hops)} = {accuracy:.1f}%")
                
                if accuracy == 100:
                    print("✓ PERFECT PREDICTION!")
                elif accuracy > 50:
                    print("⚠ Partial success - may need more observations")
                else:
                    print("✗ Reconstruction failed")
                    
            except Exception as e:
                print(f"✗ Error: {e}")
                import traceback
                traceback.print_exc()
        
        # Now test the full Gold Code reconstruction
        # (both LFSRs need to be reconstructed)
        print(f"\n{'='*70}")
        print("Full Gold Code Reconstruction (Both LFSRs)")
        print(f"{'='*70}")
        print("\nNote: Gold codes XOR two PN sequences.")
        print("To fully reconstruct, we'd need to observe the individual")
        print("PN sequences or know one of the polynomials.")
        print("\nWhat we CAN do: Reconstruct the composite sequence and predict hops!")

if __name__ == "__main__":
    test_on_dataset()
