import h5py
import os

def patch():
    dataset_path = os.path.join('data', 'synthetic', 'classification_dataset_stft_random_deg10_python.h5')
    
    if not os.path.exists(dataset_path):
        print(f"Dataset not found at {dataset_path}")
        return

    print(f"Patching attributes for {dataset_path}...")
    
    try:
        with h5py.File(dataset_path, 'r+') as f:
            # Add missing attributes
            # These values are constant across our project
            f.attrs['numHopsPerSignal'] = 256
            f.attrs['polynomial_degree'] = 10
            
            print("Attributes 'numHopsPerSignal' and 'polynomial_degree' added.")
            
            # Check if Set_ID exists (it should if the big generation finished)
            if 'Set_ID' in f:
                print("Confirmed 'Set_ID' dataset exists.")
            else:
                print("WARNING: 'Set_ID' dataset NOT found. Did generation finish?")
                
    except Exception as e:
        print(f"Error patching: {e}")

if __name__ == "__main__":
    patch()
