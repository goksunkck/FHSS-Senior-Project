
# --- train_frequency_predictor.py ---
#
# Trains a CNN-LSTM to predict the next hop frequency index using PyTorch.
# Optimized for memory efficiency (Lazy Loading) and performance.
#
# Requires:
#   - data/synthetic/prepared_prediction_sequences.mat (index + label sets)
#   - data/synthetic/classification_dataset_stft_awgn.mat (STFT tensors)

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
## TODO: Make sure to find the minimum and maximum values of the STFT data
##       and use them to normalize the data.
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

class SequenceDataset(Dataset):
    """
    Lazy-loading Dataset for large STFT data.
    Does NOT load everything into RAM. Reads from HDF5 on demand.
    Supports both MATLAB v7.3 (Object References) and Python HDF5 (Direct Arrays).
    """
    def __init__(self, h5_filepath, sequence_starts, labels, lookback_window, norm_min, norm_range, cache_size=0):
        self.h5_filepath = h5_filepath
        self.sequence_starts = sequence_starts
        self.labels = labels
        self.lookback_window = lookback_window
        self.norm_min = norm_min
        self.norm_range = norm_range
        self.pdf = None # File handle, initialized in __getitem__ to be picklable
        self.mode = 'unknown'
        
    def __len__(self):
        return len(self.sequence_starts)

    def __getitem__(self, idx):
        if self.pdf is None:
             self.pdf = h5py.File(self.h5_filepath, 'r')
             self.x_refs = self.pdf['X_data'] # Get reference to the dataset
             
             # Check Mode: MATLAB Cell Array (Refs) vs Python Tensor (Direct)
             if self.x_refs.dtype == 'object' or hasattr(self.x_refs.dtype, 'ref_dtype'):
                 self.mode = 'matlab_refs'
                 # Check if transposed (1, N)
                 if self.x_refs.shape[0] == 1 and self.x_refs.shape[1] > 1:
                     self.is_transposed = True
                 else:
                     self.is_transposed = False
             else:
                 self.mode = 'python_direct'
                 # Python shape: (N, 1024, 3)

        start_idx_matlab = self.sequence_starts[idx] # Metadata is 0-based from Python prep, 1-based from MATLAB
        label = self.labels[idx]
        
        # Adjust index based on source
        # If metadata came from Python prepare_sequences.py, it's already 0-based.
        # If came from MATLAB, it's 1-based.
        # We can detect this via filename or an attribute.
        # However, `sequence_starts` usually passed here is whatever was loaded.
        # Let's assume the LOADER handles the offset, or we check if we are reading Python file.
        
        # Heuristic: If we are in 'python_direct' mode, we assume 0-based indexing from prepare script.
        # But wait, `train()` function loads `prepared_prediction_sequences_python.h5`.
        # `prepare_sequences.py` saved `sequence_starts` as `start_idx + j` which is 0-based relative to X_data.
        # So for Python mode, start_idx is CORRECT.
        # For MATLAB mode, we subtracted 1 in the old code.
        
        start_idx = int(start_idx_matlab)
        if self.mode == 'matlab_refs':
             start_idx = start_idx - 1 # 1-based to 0-based
        
        # Safety check
        if start_idx < 0: start_idx = 0
             
        frames = []
        
        # Batch read optimization for Python mode?
        # Reading slice [start : start+LB] is faster than loop.
        if self.mode == 'python_direct':
            # shape (LB, 1024, 3)
            # x_refs is the dataset itself
            chunk = self.x_refs[start_idx : start_idx + self.lookback_window]
            
            # Normalize
            chunk = (chunk - self.norm_min) / self.norm_range
            
            # (LB, 1024, 3) -> Need (LB, C, H, W) -> (LB, 1, 1024, 3)
            # PyTorch expects C before H,W.
            # Currently (Freq, Time)? 1024x3.
            # We treat Freq=H, Time=W.
            # So (LB, H, W). Add C=1.
            # (LB, 1, 1024, 3)
            
            # Convert to Tensor
            # Optimization: Do it all at once
            chunk_tensor = torch.from_numpy(chunk).float()
            # (LB, 1024, 3) -> (LB, 3, 1024)? 
            # Train script expected (LB, 1, 1024, 3) or (LB, C, H, W).
            # Previous code:
            # img: (1024, 3, 1) -> permute(2,0,1) -> (1, 1024, 3)
            # So channel is dim 0 (of 3).
            # My Python data is (1024, 3).
            # So I want (1, 1024, 3).
            
            chunk_tensor = chunk_tensor.unsqueeze(1) # (LB, 1, 1024, 3)
            
            sequence = chunk_tensor
            
        else:
            seq_range = range(start_idx, start_idx + self.lookback_window)
            for data_idx in seq_range:
                # dereference object ref
                if self.is_transposed:
                    ref = self.x_refs[0, data_idx]
                else:
                    ref = self.x_refs[data_idx][0]
                
                img = self.pdf[ref][()] # Load specific frame
                
                # Format: Matlab usually saves as (Freq, Time). 
                # We want (C, H, W) = (1, 1024, 3) 
                if img.ndim == 2:
                    img = np.expand_dims(img, axis=-1) # (1024, 3, 1) assuming Freq x Time x C
                
                # Transpose if needed to match (H, W, C)
                if img.shape[0] == 3 and img.shape[1] > 3:
                    img = np.transpose(img, (1, 0, 2))
                    
                # Normalize
                img = (img - self.norm_min) / self.norm_range
                
                # To Tensor (H, W, C) -> (C, H, W)
                img_tensor = torch.from_numpy(img).permute(2, 0, 1).float()
                frames.append(img_tensor)

            sequence = torch.stack(frames) # (Seq, C, H, W)
            
        return sequence, torch.tensor(label, dtype=torch.long)

def train():
    print("--- 1. Load prepared metadata ---")
    # Check for Python version first
    prep_file_python = 'data/synthetic/prepared_prediction_sequences_python.h5'
    prep_file_mat = 'data/synthetic/prepared_prediction_sequences.mat'
    
    prep_file = prep_file_python if os.path.exists(prep_file_python) else prep_file_mat
    print(f"Loading metadata from: {prep_file}")
    
    if not os.path.exists(prep_file):
        print(f"Error: {prep_file} not found.")
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
    
    # --- 2. Load Normalization Stats ---
    # Try to load from the dataset file itself, fallback to fixed conservative values
    g_min = -50.0
    g_max = 50.0
    
    try:
        # We need to open the dataset file to check for global stats
        # h5py is preferred for -v7.3 mat files
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

    # --- 4. Setup Datasets ---
    # Input size check (load one frame)
    input_size = (1024, 3, 1) # Assumed default
    
    batch_size = 64 # Reverted to 64 (Faster IO interleaving)
    
    # We pass the FILENAME, not the data
    train_dataset = SequenceDataset(dataset_filename, sequence_starts_train, YTrain_mapped, lookback_window, g_min, g_range)
    val_dataset = SequenceDataset(dataset_filename, sequence_starts_validation, YValidation_mapped, lookback_window, g_min, g_range)

    # Num_workers > 0 for parallel loading
    # Windows requires if __name__ == '__main__' check which we have.
    # However, testing showed num_workers=4 is SLOWER on Windows/HDF5 due to overhead.
    # Reverting to 0 (Main Process) but keeping high Batch Size.
    train_loader = DataLoader(train_dataset, batch_size=batch_size, shuffle=True, num_workers=0, pin_memory=True)
    val_loader = DataLoader(val_dataset, batch_size=batch_size, shuffle=False, num_workers=0, pin_memory=True)

    # --- 5. Model ---
    class CNNLSTM(nn.Module):
        def __init__(self, input_size, num_classes):
            super(CNNLSTM, self).__init__()
            # Input: (batch, seq, C, H, W) -> (64, 15, 1, 1024, 3)
            # CNN expects (N, C, H, W)
            
            self.cnn = nn.Sequential(
                nn.Conv2d(1, 16, kernel_size=3, padding='same'),
                nn.BatchNorm2d(16),
                nn.ReLU(),
                nn.MaxPool2d(kernel_size=(2, 1), stride=(2, 1)), # (16, 512, 3)

                nn.Conv2d(16, 32, kernel_size=3, padding='same'),
                nn.BatchNorm2d(32),
                nn.ReLU(),
                nn.MaxPool2d(kernel_size=(2, 1), stride=(2, 1)), # (32, 256, 3)
                
                # CHANGED: AdaptiveMaxPool for better feature retention
                nn.AdaptiveMaxPool2d((10, 1)), # (32, 10, 1)
                nn.Flatten() # 320 features
            )
            
            # LSTM input size
            cnn_out_size = 32 * 10 * 1 
            
            # UPGRADED: 2 Layers, 256 Hidden Units
            # "num_layers=2" allows it to learn more complex hierarchical temporal features
            self.lstm = nn.LSTM(
                input_size=cnn_out_size, 
                hidden_size=256, 
                num_layers=2, 
                batch_first=True,
                dropout=0.5 # Dropout between LSTM layers
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
    
    model = CNNLSTM(input_size, num_classes).to(device)
    
    # --- 6. Train ---
    epochs = 20
    criterion = nn.CrossEntropyLoss()
    # Changed to AdamW with weight_decay and explicit betas
    optimizer = torch.optim.AdamW(model.parameters(), lr=1e-4, betas=(0.9, 0.999), weight_decay=1e-4) # "Thumb rule" betas
    
    # Scheduler: Reduce LR if val_loss doesn't improve for 2 epochs
    scheduler = torch.optim.lr_scheduler.ReduceLROnPlateau(optimizer, mode='min', factor=0.5, patience=2)
    
    history = {'train_loss': [], 'val_loss': [], 'val_accuracy': []}
    
    print("Starting training with Upgraded Model (2-Layer LSTM, 256 Hidden)...")
    
    import time
    
    best_val_loss = float('inf')
    
    import time
    
    for epoch in range(epochs):
        model.train()
        running_loss = 0.0
        
        epoch_start_time = time.time()
        num_batches = len(train_loader)
        
        for i, (sequences, labels) in enumerate(train_loader):
            step_start = time.time()
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
        
        # Step Scheduler
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
    
    metadata = {
        'normalization': {'min': float(g_min), 'range': float(g_range)},
        'history': history,
        'class_values': class_values.tolist()
    }
    with open('models/model_metadata.json', 'w') as f:
        json.dump(metadata, f)
        
    print("Done.")

if __name__ == '__main__':
    train()
