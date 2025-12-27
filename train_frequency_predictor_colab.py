
# --- train_frequency_predictor_colab.py ---
#
# Optimized for Google Colab / High-RAM Environments.
# LOADS ENTIRE DATASET INTO RAM for maximum training speed.
#
# Requires ~12GB+ RAM depending on dataset size.

import json
import numpy as np
import scipy.io
import torch
import torch.nn as nn
from torch.utils.data import Dataset, DataLoader
import h5py
import os
import matplotlib.pyplot as plt
import random

def load_mat_file(filepath, variable_name):
    """Loads a specific variable from a .mat file."""
    try:
        return scipy.io.loadmat(filepath)
    except NotImplementedError:
        import h5py
        with h5py.File(filepath, 'r') as f:
            data = {}
            for k, v in f.items():
                if isinstance(v, h5py.Dataset):
                    value = v[()]
                    if isinstance(value, np.ndarray):
                        if value.shape == (1, 1):
                            data[k] = value[0, 0]
                        else:
                            data[k] = value
                    else:
                        data[k] = value
            return data

class RAMDataset(Dataset):
    """
    High-Performance Dataset that loads EVERYTHING into RAM at init.
    Ideal for Colab.
    """
    def __init__(self, h5_filepath, sequence_starts, labels, lookback_window, norm_min, norm_range):
        print(f"preload: Loading dataset from {h5_filepath} into RAM...")
        self.sequence_starts = sequence_starts
        self.labels = labels
        self.lookback_window = lookback_window
        
        # Load the entire HDF5 file content into a CPU Tensor or Numpy Array
        with h5py.File(h5_filepath, 'r') as f:
            # Check if transposed
            x_refs = f['X_data']
            
            # We iterate and load all referenced frames now.
            # Warning: This loop might take a minute, but training will be blazing fast.
            # If the HDF5 structure is simple (just a big array), we load it directly. 
            # If it uses cell object references (MATLAB -v7.3), we must dereference.
            
            # --- Optimized Bulk Loading ---
            # Since checking 60,000 references is slow in Python, we try to load raw
            # If X_data is a cell array in MATLAB, it's a Dataset of Object References in HDF5.
            
            self.frames = []
            
            # Pre-load Loop
            # We only load frames that are actually USED in the sequences to save some RAM
            # But the user asked for "Max Potential", so let's load what we need.
            
            # Flatten indices needed
            needed_indices = set()
            for start in sequence_starts:
                # Matlab 1-based -> Python 0-based
                s_idx = start - 1
                for i in range(lookback_window):
                    needed_indices.add(s_idx + i)
            
            # Convert to sorted list for sequential access (faster HDD seek)
            needed_indices = sorted(list(needed_indices))
            
            # Create a localized cache/map: global_index -> tensor
            self.cache = {}
            
            # Determine dereferencing method once
            is_transposed = (x_refs.shape[0] == 1 and x_refs.shape[1] > 1)
            
            for idx in needed_indices:
                if is_transposed:
                    ref = x_refs[0, idx]
                else:
                    ref = x_refs[idx][0]
                
                img = f[ref][()]
                
                # Reshape/Transposes (Same logic as lazy loader)
                if img.ndim == 2:
                    img = np.expand_dims(img, axis=-1)
                
                if img.shape[0] == 3 and img.shape[1] > 3:
                     # (3, 1024, 1) -> (1024, 3, 1)
                    img = np.transpose(img, (1, 0, 2))
                    
                # Normalize immediately to save GPU compute
                img = (img - norm_min) / norm_range
                
                # Store as Half precision to save RAM (optional, but good for Colab)
                # keeping float32 for safety
                self.cache[idx] = torch.from_numpy(img).permute(2, 0, 1).float()
                
        print(f"preload: Loaded {len(self.cache)} frames into RAM.")

    def __len__(self):
        return len(self.sequence_starts)

    def __getitem__(self, idx):
        start_idx_matlab = self.sequence_starts[idx]
        label = self.labels[idx]
        start_idx_python = start_idx_matlab - 1
        
        frames = []
        for i in range(self.lookback_window):
            # Fetch from RAM cache
            frames.append(self.cache[start_idx_python + i])
            
        sequence = torch.stack(frames)
        return sequence, torch.tensor(label, dtype=torch.long)

def train():
    print("--- 1. Load prepared metadata ---")
    prep_file = 'data/synthetic/prepared_prediction_sequences.mat'
    
    # Check file existence
    if not os.path.exists(prep_file):
         print(f"Error: {prep_file} not found. Please upload data.")
         return

    prep_data = load_mat_file(prep_file, 'prep_data')
    sequence_starts_train = prep_data['sequence_starts_train'].flatten()
    sequence_starts_validation = prep_data['sequence_starts_validation'].flatten()
    YTrain = prep_data['YTrain'].flatten()
    YValidation = prep_data['YValidation'].flatten()
    lookback_window = int(prep_data['lookback_window'])
    class_values = prep_data['class_values'].flatten()
    
    dataset_filename = prep_data['dataset_filename']
    if isinstance(dataset_filename, np.ndarray):
        dataset_filename = "".join(map(chr, dataset_filename.flatten()))

    print(f"Dataset path: {dataset_filename}")
    
    # --- 2. Normalization ---
    # Try to load from the dataset file itself, fallback to fixed conservative values
    g_min = -50.0
    g_max = 50.0
    
    try:
        with h5py.File(dataset_filename, 'r') as f:
            if 'global_min' in f and 'global_max' in f:
                g_min = float(f['global_min'][0, 0])
                g_max = float(f['global_max'][0, 0])
                print(f"Loaded normalization stats from file: Min={g_min:.2f}, Max={g_max:.2f}")
            else:
                 print("Global stats not found in file. Using default conservative values.")
    except Exception as e:
        print(f"Warning: Could not read normalization stats ({e}). Using defaults.")

    print(f"Using normalization: Min={g_min}, Max={g_max}")
    g_range = g_max - g_min

    # --- 3. Process Labels ---
    class_map = {val: i for i, val in enumerate(class_values)}
    YTrain_mapped = np.array([class_map[y] for y in YTrain], dtype=np.int64)
    YValidation_mapped = np.array([class_map[y] for y in YValidation], dtype=np.int64)
    num_classes = len(class_values)

    # --- 4. Setup Datasets (RAM LOAD) ---
    print("Initializing RAM Datasets (This may take 1-2 mins)...")
    
    # Use the RAMDataset class
    train_dataset = RAMDataset(dataset_filename, sequence_starts_train, YTrain_mapped, lookback_window, g_min, g_range)
    val_dataset = RAMDataset(dataset_filename, sequence_starts_validation, YValidation_mapped, lookback_window, g_min, g_range)

    # High Performance Loader Settings
    # Colab has 2-4 cores. num_workers=4 is usually safe.
    batch_size = 256 # High batch size for GPU saturation
    
    train_loader = DataLoader(train_dataset, batch_size=batch_size, shuffle=True, num_workers=4, pin_memory=True)
    val_loader = DataLoader(val_dataset, batch_size=batch_size, shuffle=False, num_workers=4, pin_memory=True)

    # --- 5. Model (Same upgraded architecture) ---
    class CNNLSTM(nn.Module):
        def __init__(self, input_size, num_classes):
            super(CNNLSTM, self).__init__()
            # CNN
            self.cnn = nn.Sequential(
                nn.Conv2d(1, 16, kernel_size=3, padding='same'),
                nn.BatchNorm2d(16),
                nn.ReLU(),
                nn.MaxPool2d(kernel_size=(2, 1), stride=(2, 1)),

                nn.Conv2d(16, 32, kernel_size=3, padding='same'),
                nn.BatchNorm2d(32),
                nn.ReLU(),
                nn.MaxPool2d(kernel_size=(2, 1), stride=(2, 1)),
                
                nn.AdaptiveMaxPool2d((10, 1)), 
                nn.Flatten() 
            )
            
            cnn_out_size = 32 * 10 * 1 
            
            # LSTM (2-Layer, 256 Hidden)
            self.lstm = nn.LSTM(
                input_size=cnn_out_size, 
                hidden_size=256, 
                num_layers=2, 
                batch_first=True,
                dropout=0.5
            )
            
            self.dropout = nn.Dropout(0.5)
            self.fc = nn.Linear(256, num_classes)

        def forward(self, x):
            batch_size, seq_len, C, H, W = x.size()
            c_in = x.view(batch_size * seq_len, C, H, W)
            c_out = self.cnn(c_in)
            r_in = c_out.view(batch_size, seq_len, -1)
            _, (h_n, _) = self.lstm(r_in)
            out = self.dropout(h_n[-1])
            out = self.fc(out)
            return out

    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    print(f"Using device: {device}")
    
    # Determine Input Size from cache
    # Get one random frame
    first_key = next(iter(train_dataset.cache))
    sample_frame = train_dataset.cache[first_key] # (C, H, W)
    input_size = sample_frame.shape # e.g. (1, 1024, 3)
    
    model = CNNLSTM(input_size, num_classes).to(device)
    
    # --- 6. Train ---
    epochs = 50 # Increased for Colab
    criterion = nn.CrossEntropyLoss()
    
    # AdamW + Scheduler
    optimizer = torch.optim.AdamW(model.parameters(), lr=1e-4, betas=(0.9, 0.999), weight_decay=1e-4) # "Thumb rule" betas
    scheduler = torch.optim.lr_scheduler.ReduceLROnPlateau(optimizer, mode='min', factor=0.5, patience=3)
    
    history = {'train_loss': [], 'val_loss': [], 'val_accuracy': []}
    
    print("Starting High-Performance Training...")
    
    import time
    best_val_loss = float('inf')
    
    for epoch in range(epochs):
        model.train()
        running_loss = 0.0
        
        epoch_start_time = time.time()
        num_batches = len(train_loader)
        
        for i, (sequences, labels) in enumerate(train_loader):
            sequences = sequences.to(device)
            labels = labels.to(device)
            
            optimizer.zero_grad()
            outputs = model(sequences)
            loss = criterion(outputs, labels)
            loss.backward()
            optimizer.step()
            
            running_loss += loss.item()
            
            # Print status every 10 batches
            if (i + 1) % 10 == 0 or (i + 1) == num_batches:
                elapsed = time.time() - epoch_start_time
                batches_done = i + 1
                avg_time_per_batch = elapsed / batches_done
                remaining_batches = num_batches - batches_done
                eta_seconds = avg_time_per_batch * remaining_batches
                
                print(f"Epoch [{epoch+1}/{epochs}] Step [{batches_done}/{num_batches}] "
                      f"Loss: {loss.item():.4f} "
                      f"| Elapsed: {elapsed:.0f}s | ETA: {eta_seconds:.0f}s")
            
        avg_train_loss = running_loss / len(train_loader)
        
        # Validation
        model.eval()
        correct = 0
        total = 0
        val_loss = 0
        with torch.no_grad():
            for sequences, labels in val_loader:
                sequences = sequences.to(device)
                labels = labels.to(device)
                outputs = model(sequences)
                val_loss += criterion(outputs, labels).item()
                _, predicted = torch.max(outputs.data, 1)
                total += labels.size(0)
                correct += (predicted == labels).sum().item()
        
        avg_val_loss = val_loss / len(val_loader)
        acc = 100 * correct / total
        
        scheduler.step(avg_val_loss)
        
        print(f"Epoch [{epoch+1}/{epochs}] Train Loss: {avg_train_loss:.4f}, Val Loss: {avg_val_loss:.4f}, Acc: {acc:.2f}%")
        
        history['train_loss'].append(avg_train_loss)
        history['val_loss'].append(avg_val_loss)
        history['val_accuracy'].append(acc)
        
        # Save Best Model
        if avg_val_loss < best_val_loss:
            best_val_loss = avg_val_loss
            print(f"  -> Validation Loss Improved. Saving model...")
            if not os.path.exists('models'): os.makedirs('models')
            torch.save(model.state_dict(), 'models/trained_frequency_predictor.pth')
            
    print("Done.")

if __name__ == '__main__':
    train()
