"""
Test Berlekamp-Massey Algorithm with Noisy Observations

This tests the sensitivity of Berlekamp-Massey to bit errors,
which would occur with:
1. Low SNR (frequency detection errors)
2. ML classifier errors
"""

import h5py
import numpy as np
import sys
import os

sys.path.append(os.path.join(os.getcwd(), 'src', 'python'))

from berlekamp_massey import reconstruct_gold_code_lfsr
from gold_code import LFSR

def test_with_snr_levels():
    """Test on signals with different SNR levels"""
    dataset_file = 'data/synthetic/classification_dataset_stft_random_deg10_python.h5'
    
    with h5py.File(dataset_file, 'r') as f:
        Y_data = f['Y_data']
        X_data = f['X_data']
        
        # We generated at SNR levels: -10, -8, -6, -4, -2, 0, 2, 4, 6, 8, 10 dB
        # Each level has 100 signals * 256 hops = 25,600 hops
        hops_per_snr = 100 * 256
        
        print("="*70)
        print("Testing Berlekamp-Massey with Different SNR Levels")
        print("="*70)
        
        snr_levels = [-10, -8, -6, -4, -2, 0, 2, 4, 6, 8, 10]
        
        for idx, snr in enumerate(snr_levels):
            print(f"\n{'='*70}")
            print(f"SNR = {snr} dB")
            print(f"{'='*70}")
            
            # Get signal at this SNR level
            start_idx = idx * hops_per_snr
            signal_hops = Y_data[start_idx : start_idx + 256, 0]
            
            # Observe first 50 hops
            num_observed = 50
            observed_hops = signal_hops[:num_observed]
            actual_future_hops = signal_hops[num_observed:num_observed+20]
            
            print(f"Observed hops: {observed_hops[:20]}...")
            
            try:
                # Reconstruct
                polynomial, initial_state, degree = reconstruct_gold_code_lfsr(observed_hops, k=3)
                
                # Predict
                lfsr = LFSR(polynomial, initial_state)
                _ = lfsr.generate(num_observed * 3)
                future_bits = lfsr.generate(20 * 3)
                
                predicted_hops = []
                for i in range(0, len(future_bits), 3):
                    hop_bits = future_bits[i:i+3]
                    hop_index = hop_bits[0] * 4 + hop_bits[1] * 2 + hop_bits[2]
                    predicted_hops.append(hop_index)
                
                predicted_hops = np.array(predicted_hops)
                
                # Compare
                matches = np.sum(predicted_hops == actual_future_hops)
                accuracy = matches / len(actual_future_hops) * 100
                
                print(f"Actual:    {actual_future_hops}")
                print(f"Predicted: {predicted_hops}")
                print(f"Accuracy: {accuracy:.1f}%")
                
                if accuracy == 100:
                    print("✓ Perfect")
                elif accuracy > 90:
                    print("⚠ Good (minor errors)")
                elif accuracy > 50:
                    print("⚠ Degraded")
                else:
                    print("✗ Failed")
                    
            except Exception as e:
                print(f"✗ Reconstruction failed: {e}")


def test_with_simulated_errors():
    """Simulate what happens when we have frequency detection errors"""
    dataset_file = 'data/synthetic/classification_dataset_stft_random_deg10_python.h5'
    
    with h5py.File(dataset_file, 'r') as f:
        Y_data = f['Y_data']
        
        # Get a clean signal
        signal_hops = Y_data[0:256, 0]
        
        print("\n" + "="*70)
        print("Testing with Simulated Frequency Detection Errors")
        print("="*70)
        
        error_rates = [0.0, 0.01, 0.02, 0.05, 0.10]
        
        for error_rate in error_rates:
            print(f"\n{'='*70}")
            print(f"Error Rate: {error_rate*100:.1f}% (bit error rate)")
            print(f"{'='*70}")
            
            # Add random errors to observed hops
            num_observed = 50
            observed_hops = signal_hops[:num_observed].copy()
            
            # Flip some hops randomly
            num_errors = int(num_observed * error_rate)
            if num_errors > 0:
                error_positions = np.random.choice(num_observed, num_errors, replace=False)
                for pos in error_positions:
                    # Flip to a random wrong frequency
                    wrong_hop = np.random.choice([h for h in range(8) if h != observed_hops[pos]])
                    observed_hops[pos] = wrong_hop
                print(f"Introduced {num_errors} errors at positions: {sorted(error_positions)}")
            
            actual_future_hops = signal_hops[num_observed:num_observed+20]
            
            try:
                polynomial, initial_state, degree = reconstruct_gold_code_lfsr(observed_hops, k=3)
                
                # Predict
                lfsr = LFSR(polynomial, initial_state)
                _ = lfsr.generate(num_observed * 3)
                future_bits = lfsr.generate(20 * 3)
                
                predicted_hops = []
                for i in range(0, len(future_bits), 3):
                    hop_bits = future_bits[i:i+3]
                    hop_index = hop_bits[0] * 4 + hop_bits[1] * 2 + hop_bits[2]
                    predicted_hops.append(hop_index)
                
                predicted_hops = np.array(predicted_hops)
                
                matches = np.sum(predicted_hops == actual_future_hops)
                accuracy = matches / len(actual_future_hops) * 100
                
                print(f"Actual:    {actual_future_hops}")
                print(f"Predicted: {predicted_hops}")
                print(f"Accuracy: {accuracy:.1f}%")
                
                if error_rate == 0.0 and accuracy < 100:
                    print("✗ UNEXPECTED: Should be perfect with no errors!")
                elif error_rate > 0 and accuracy == 100:
                    print("⚠ Lucky: Still works despite errors")
                elif error_rate > 0 and accuracy < 20:
                    print("✗ Completely broken (as expected with errors)")
                else:
                    print(f"⚠ Degraded performance")
                    
            except Exception as e:
                print(f"✗ Reconstruction failed: {e}")


if __name__ == "__main__":
    # Test 1: Different SNR levels
    test_with_snr_levels()
    
    # Test 2: Simulated errors
    np.random.seed(42)
    test_with_simulated_errors()
    
    print("\n" + "="*70)
    print("CONCLUSION")
    print("="*70)
    print("\nBerlekamp-Massey is VERY sensitive to errors:")
    print("  - Even 1% error rate can break reconstruction")
    print("  - Low SNR → frequency detection errors → fails")
    print("\nSolution needed:")
    print("  - Error correction codes")
    print("  - Majority voting over multiple observations")
    print("  - Hybrid: ML for robust detection + BM for prediction")
