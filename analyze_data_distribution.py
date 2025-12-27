import h5py
import numpy as np
import matplotlib.pyplot as plt

def analyze():
    dataset_path = 'data/synthetic/classification_dataset_stft_random_deg10_python.h5'
    
    with h5py.File(dataset_path, 'r') as f:
        X = f['X_data']
        # Sample random chunks
        print("Sampling data...")
        indices = np.random.choice(X.shape[0], 1000, replace=False)
        indices.sort() # h5py requires sorted indices
        samples = X[indices] # (1000, 1024, 3)
        
        flat = samples.flatten()
        
        print(f"Global Stats (File): Min={f.attrs['global_min'] if 'global_min' in f.attrs else 'N/A'}, Max={f.attrs['global_max'] if 'global_max' in f.attrs else 'N/A'}")
        print(f"Sample Stats: Min={flat.min()}, Max={flat.max()}, Mean={flat.mean()}, Std={flat.std()}")
        print(f"Percentiles: 1%={np.percentile(flat, 1)}, 5%={np.percentile(flat, 5)}, 50%={np.median(flat)}, 95%={np.percentile(flat, 95)}, 99%={np.percentile(flat, 99)}")
        
        plt.figure()
        plt.hist(flat, bins=100)
        plt.title("Histogram of STFT Log Magnitudes (dB)")
        plt.xlabel("dB")
        plt.ylabel("Count")
        plt.savefig('results/data_histogram.png')
        print("Saved histogram to results/data_histogram.png")

if __name__ == "__main__":
    analyze()
