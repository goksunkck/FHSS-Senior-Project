
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
        self.norm_min = norm_min
        self.norm_range = norm_range
        
        # Load the entire HDF5 file content
        with h5py.File(h5_filepath, 'r') as f:
            x_data = f['X_data']
            
            # Check Mode: MATLAB Cell Array (Refs) vs Python Tensor (Direct)
            if x_data.dtype == 'object' or hasattr(x_data.dtype, 'ref_dtype'):
                 self.mode = 'matlab_refs'
                 print("Mode: Legacy MATLAB (Object References)")
            else:
                 self.mode = 'python_direct'
                 print("Mode: Python Direct (HDF5 Array)")

            if self.mode == 'python_direct':
                # --- FAST PATH ---
                # Load everything directly
                # Shape: (TotalHops, 1024, 3)
                print("Loading raw data into RAM...")
                raw_data = x_data[:] 
                
                # Normalize
                raw_data = (raw_data - norm_min) / norm_range
                
                # Convert to Tensor (N, 1024, 3)
                self.data_tensor = torch.from_numpy(raw_data).float()
                
                # Reshape to (N, 1, 1024, 3) for CNN
                self.data_tensor = self.data_tensor.unsqueeze(1)
                
                print(f"preload: Loaded {self.data_tensor.shape} tensor into RAM.")
                
            else:
                # --- LEGACY PATH (Refs) ---
                x_refs = x_data
                self.frames = []
                self.cache = {}
                
                needed_indices = set()
                for start in sequence_starts:
                    # Matlab 1-based -> Python 0-based
                    s_idx = start - 1
                    for i in range(lookback_window):
                        needed_indices.add(s_idx + i)
                
                needed_indices = sorted(list(needed_indices))
                is_transposed = (x_refs.shape[0] == 1 and x_refs.shape[1] > 1)
                
                for idx in needed_indices:
                    if is_transposed:
                        ref = x_refs[0, idx]
                    else:
                        ref = x_refs[idx][0]
                    
                    img = f[ref][()]
                    
                    if img.ndim == 2:
                        img = np.expand_dims(img, axis=-1)
                    
                    if img.shape[0] == 3 and img.shape[1] > 3:
                        img = np.transpose(img, (1, 0, 2))
                        
                    img = (img - norm_min) / norm_range
                    self.cache[idx] = torch.from_numpy(img).permute(2, 0, 1).float()
                    
                print(f"preload: Loaded {len(self.cache)} frames into RAM.")

    def __len__(self):
        return len(self.sequence_starts)

    def __getitem__(self, idx):
        # Adjust start index based on mode
        # MATLAB metadata is 1-based. Python metadata is 0-based.
        # BUT: the 'sequence_starts' passed here comes from loaded metadata.
        # We need to rely on the loader logic to give us the raw value.
        # If Python mode, we assume 0-based. If MATLAB, 1-based.
        
        # NOTE: Logic below assumes we subtract 1 for MATLAB. 
        # For Python: 'prepare_sequences.py' saved 0-based indices.
        
        start_raw = self.sequence_starts[idx]
        label = self.labels[idx]
        
        if self.mode == 'python_direct':
            # Index is 0-based, direct slice
            start = int(start_raw)
            # Slice range
            end = start + self.lookback_window
            
            # Direct tensor slice (Fast)
            sequence = self.data_tensor[start:end]
            
            return sequence, torch.tensor(label, dtype=torch.long)
            
        else:
            # Legacy
            start_idx_python = start_raw - 1
            frames = []
            for i in range(self.lookback_window):
                frames.append(self.cache[start_idx_python + i])
            sequence = torch.stack(frames)
            return sequence, torch.tensor(label, dtype=torch.long)

def train():
    print("--- 1. Load prepared metadata ---")
    
    # Priority: Python H5 > MATLAB MAT
    prep_file_python = 'data/synthetic/prepared_prediction_sequences_python.h5'
    prep_file_legacy = 'data/synthetic/prepared_prediction_sequences.mat'
    
    prep_data = {}
    is_python_source = False
    
    if os.path.exists(prep_file_python):
        print(f"Loading metadata from: {prep_file_python}")
        with h5py.File(prep_file_python, 'r') as f:
            prep_data['sequence_starts_train'] = f['sequence_starts_train'][:]
            prep_data['sequence_starts_validation'] = f['sequence_starts_validation'][:]
            prep_data['YTrain'] = f['YTrain'][:]
            prep_data['YValidation'] = f['YValidation'][:]
            # Scalar attributes
            if 'lookback_window' in f.attrs:
                prep_data['lookback_window'] = int(f.attrs['lookback_window'])
            elif 'lookback_window' in f:
                prep_data['lookback_window'] = int(f['lookback_window'][0])
            else:
                raise KeyError("lookback_window not found in metadata")

            # String handling
            if 'dataset_filename' in f.attrs:
                dset_name = f.attrs['dataset_filename']
                if isinstance(dset_name, (bytes, np.bytes_)):
                     prep_data['dataset_filename'] = dset_name.decode('utf-8')
                else:
                     prep_data['dataset_filename'] = str(dset_name)
            elif 'dataset_filename' in f:
                dset_name = f['dataset_filename'][:]
                if isinstance(dset_name[0], (bytes, np.bytes_)):
                     prep_data['dataset_filename'] = dset_name[0].decode('utf-8')
                else:
                     prep_data['dataset_filename'] = "".join([chr(c) for c in dset_name.flatten()])
            else:
                # Fallback or error
                # If we are in python mode, we know the filename usually
                prep_data['dataset_filename'] = 'data/synthetic/classification_dataset_stft_random_deg10_python.h5'
            
            # Class values
            prep_data['class_values'] = f['class_values'][:]
            is_python_source = True
            
    elif os.path.exists(prep_file_legacy):
        print(f"Loading metadata from: {prep_file_legacy}")
        data = load_mat_file(prep_file_legacy, 'prep_data')
        prep_data['sequence_starts_train'] = data['sequence_starts_train'].flatten()
        prep_data['sequence_starts_validation'] = data['sequence_starts_validation'].flatten()
        prep_data['YTrain'] = data['YTrain'].flatten()
        prep_data['YValidation'] = data['YValidation'].flatten()
        prep_data['lookback_window'] = int(data['lookback_window'])
        prep_data['class_values'] = data['class_values'].flatten()
        
        d_name = data['dataset_filename']
        if isinstance(d_name, np.ndarray):
            prep_data['dataset_filename'] = "".join(map(chr, d_name.flatten()))
        else:
            prep_data['dataset_filename'] = str(d_name)
    else:
        print("Error: No prepared metadata found.")
        return

    sequence_starts_train = prep_data['sequence_starts_train']
    sequence_starts_validation = prep_data['sequence_starts_validation']
    YTrain = prep_data['YTrain']
    YValidation = prep_data['YValidation']
    lookback_window = prep_data['lookback_window']
    class_values = prep_data['class_values']
    dataset_filename = prep_data['dataset_filename']

    print(f"Dataset path: {dataset_filename}")
    
    # --- 2. Normalization ---
    # Load from file (attributes or datasets)
    g_min = -100.0
    g_max = 50.0
    
    try:
        with h5py.File(dataset_filename, 'r') as f:
            if 'global_min' in f.attrs:
                g_min = float(f.attrs['global_min'])
                g_max = float(f.attrs['global_max'])
                print(f"Loaded normalization from attributes: Min={g_min:.2f}, Max={g_max:.2f}")
            elif 'global_min' in f:
                g_min = float(f['global_min'][0])
                g_max = float(f['global_max'][0])
                print(f"Loaded normalization from datasets: Min={g_min:.2f}, Max={g_max:.2f}")
            else:
                print("Warning: Global stats not in file. Using defaults.")
    except Exception as e:
        print(f"Warning: Could not read normalization ({e}). Using defaults.")
    
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
    # Windows: num_workers=0 to prevent RAM explosion with large Dataset objects
    # Batch size reduced to fit 6GB VRAM
    batch_size = 32
    
    train_loader = DataLoader(train_dataset, batch_size=batch_size, shuffle=True, num_workers=0, pin_memory=True)
    val_loader = DataLoader(val_dataset, batch_size=batch_size, shuffle=False, num_workers=0, pin_memory=True)

    # --- 5. Model (Improved Architecture) ---
    class CNNLSTM(nn.Module):
        def __init__(self, input_size, num_classes):
            super(CNNLSTM, self).__init__()
            # Input: (1, 1024, 3)
            
            self.cnn = nn.Sequential(
                # Block 1
                nn.Conv2d(1, 32, kernel_size=3, padding='same'),
                nn.BatchNorm2d(32),
                nn.ReLU(),
                nn.MaxPool2d(kernel_size=(2, 1)), # -> (32, 512, 3)
                
                # Block 2
                nn.Conv2d(32, 64, kernel_size=3, padding='same'),
                nn.BatchNorm2d(64),
                nn.ReLU(),
                nn.MaxPool2d(kernel_size=(2, 1)), # -> (64, 256, 3)
                
                # Block 3
                nn.Conv2d(64, 128, kernel_size=3, padding='same'),
                nn.BatchNorm2d(128),
                nn.ReLU(),
                nn.MaxPool2d(kernel_size=(2, 1)), # -> (128, 128, 3)
                
                # Block 4
                nn.Conv2d(128, 128, kernel_size=3, padding='same'),
                nn.BatchNorm2d(128),
                nn.ReLU(),
                nn.MaxPool2d(kernel_size=(4, 1)), # -> (128, 32, 3)
            )
            
            # Calculate flatten size
            # Final Height = 32
            # Final Width = 3 (time dim preserved)
            # Channels = 128
            self.flatten_size = 128 * 32 * 3
            
            self.projection = nn.Linear(self.flatten_size, 512)
            self.proj_bn = nn.BatchNorm1d(512)
            self.proj_relu = nn.ReLU()
            
            # LSTM (2-Layer, 512 Hidden)
            self.lstm = nn.LSTM(
                input_size=512, 
                hidden_size=512, 
                num_layers=2, 
                batch_first=True,
                dropout=0.5
            )
            
            self.dropout = nn.Dropout(0.5)
            self.fc = nn.Linear(512, num_classes)

        def forward(self, x):
            batch_size, seq_len, C, H, W = x.size()
            c_in = x.view(batch_size * seq_len, C, H, W)
            
            c_out = self.cnn(c_in) # (B*S, 128, 32, 3)
            
            c_flat = c_out.view(batch_size * seq_len, -1)
            c_proj = self.proj_relu(self.proj_bn(self.projection(c_flat)))
            
            r_in = c_proj.view(batch_size, seq_len, -1)
            
            # LSTM
            # output: (B, S, H_out), h_n: (Layers, B, H_out)
            lstm_out, (h_n, _) = self.lstm(r_in)
            
            # Use last hidden state
            # h_n[-1] is the last layer's final state
            out = self.dropout(h_n[-1])
            out = self.fc(out)
            return out

    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    print(f"Using device: {device}")
    
    # Determine Input Size from cache or tensor
    if train_dataset.mode == 'python_direct':
        # (N, 1, 1024, 3) -> frame (1, 1024, 3)
        input_size = train_dataset.data_tensor.shape[1:] 
    else:
        first_key = next(iter(train_dataset.cache))
        sample_frame = train_dataset.cache[first_key] # (C, H, W)
        input_size = sample_frame.shape 
    
    print(f"Input Frame Size: {input_size}")
    
    model = CNNLSTM(input_size, num_classes).to(device)
    
    # --- 6. Train ---
    epochs = 50 # Increased for Colab
    criterion = nn.CrossEntropyLoss()
    
    # AdamW + Scheduler
    # Increased LR to 1e-3
    optimizer = torch.optim.AdamW(model.parameters(), lr=1e-3, betas=(0.9, 0.999), weight_decay=1e-4)
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
