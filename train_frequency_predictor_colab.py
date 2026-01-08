
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
    def __init__(self, h5_filepath, sequence_starts, labels, lookback_window):
        print(f"preload: Loading dataset from {h5_filepath} into RAM...")
        self.sequence_starts = sequence_starts
        self.labels = labels
        self.lookback_window = lookback_window
        
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
                
                # Normalize PER SAMPLE (Vectorized)
                # raw_data shape: (TotalHops, 1024, 3)
                
                print("Applying per-sample normalization for maximum contrast...")
                
                # Reshape to (N, -1) to find min/max per sample
                flat_data = raw_data.reshape(raw_data.shape[0], -1)
                
                # Get min and max per sample
                mins = flat_data.min(axis=1, keepdims=True)
                maxs = flat_data.max(axis=1, keepdims=True)
                ranges = maxs - mins
                
                # Handle zero range (prevent div by zero)
                ranges[ranges == 0] = 1.0
                
                # Expand dims back to (N, 1, 1) for broadcasting against (N, 1024, 3)
                # Actually min/max are (N, 1). We need (N, 1, 1)
                mins = mins.reshape(-1, 1, 1)
                ranges = ranges.reshape(-1, 1, 1)
                
                # Normalize
                raw_data = (raw_data - mins) / ranges
                
                print("Per-sample normalization complete.")
                
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
                        
                    # Per-sample normalization
                    i_min = img.min()
                    i_max = img.max()
                    i_range = i_max - i_min
                    if i_range == 0: i_range = 1.0
                    
                    img = (img - i_min) / i_range
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
    prep_file_python = os.path.join('data', 'synthetic', 'prepared_prediction_sequences_python.h5')
    prep_file_legacy = os.path.join('data', 'synthetic', 'prepared_prediction_sequences.mat')
    
    print(f"Checking for Python metadata at: {os.path.abspath(prep_file_python)}")
    print(f"Checking for Legacy metadata at: {os.path.abspath(prep_file_legacy)}")
    
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
    # SKIPPED: Using per-sample normalization inside RAMDataset
    print(f"Using per-sample normalization (0-1 range).")

    # --- 3. Process Labels ---
    class_map = {val: i for i, val in enumerate(class_values)}
    YTrain_mapped = np.array([class_map[y] for y in YTrain], dtype=np.int64)
    YValidation_mapped = np.array([class_map[y] for y in YValidation], dtype=np.int64)
    num_classes = len(class_values)

    # --- 4. Setup Datasets (RAM LOAD) ---
    print("Initializing RAM Datasets (This may take 1-2 mins)...")
    
    # Use the RAMDataset class
    train_dataset = RAMDataset(dataset_filename, sequence_starts_train, YTrain_mapped, lookback_window)
    val_dataset = RAMDataset(dataset_filename, sequence_starts_validation, YValidation_mapped, lookback_window)

    # High Performance Loader Settings
    # Windows: num_workers=0 to prevent RAM explosion with large Dataset objects
    # Batch size reduced to fit 6GB VRAM
    batch_size = 128
    
    train_loader = DataLoader(train_dataset, batch_size=batch_size, shuffle=True, num_workers=0, pin_memory=True)
    val_loader = DataLoader(val_dataset, batch_size=batch_size, shuffle=False, num_workers=0, pin_memory=True)

    # --- 5. Model (Improved Architecture) ---
    class CNNLSTM(nn.Module):
        def __init__(self, input_size, num_classes):
            super(CNNLSTM, self).__init__()
            
            # 1. Vision Module (Pre-trained CNN)
            self.cnn = nn.Sequential(
                nn.Conv2d(1, 16, kernel_size=(5, 3), padding=(2, 0)),
                nn.InstanceNorm2d(16), 
                nn.ReLU(),
                nn.MaxPool2d(kernel_size=(4, 1)), 
                
                nn.Conv2d(16, 32, kernel_size=(5, 1), padding='same'),
                nn.InstanceNorm2d(32), 
                nn.ReLU(),
                nn.MaxPool2d(kernel_size=(4, 1)), 
            )
            
            # CNN Output Head (to get the class ID)
            # Input: 32 * 64 * 1 = 2048
            self.cnn_fc = nn.Linear(32 * 64, num_classes) # Outputs 8 logits
            
            # 2. Symbolic Interface
            # We take the argmax of CNN (0-7) and embed it
            self.embedding = nn.Embedding(num_embeddings=num_classes, embedding_dim=32)
            
            # 3. Sequence Modeler (LSTM)
            # Input size is now just embedding_dim (32)
            self.lstm = nn.LSTM(
                input_size=32,
                hidden_size=64, # Sufficient for 20-bit state
                num_layers=2, 
                batch_first=True,
                dropout=0.2
            )
            
            self.fc = nn.Linear(64, num_classes)

        def forward(self, x):
            # x: (Batch, Seq, C, H, W)
            batch_size, seq_len, C, H, W = x.size()
            
            # Fold batch and sequence dims for CNN
            c_in = x.view(batch_size * seq_len, C, H, W)
            
            # CNN Forward (No Gradients for Vision Part)
            with torch.no_grad():
                c_out = self.cnn(c_in) # (B*S, 32, 64, 1)
                c_flat = c_out.view(c_out.size(0), -1) # (B*S, 2048)
                cnn_logits = self.cnn_fc(c_flat) # (B*S, 8)
                _, predicted_indices = torch.max(cnn_logits, dim=1) # (B*S, )
            
            # Symbolic Embedding (Gradient Flow Starts Here)
            embedded = self.embedding(predicted_indices) # (B*S, 32)
            
            # Reshape to sequence for LSTM
            lstm_input = embedded.view(batch_size, seq_len, -1) # (B, S, 32)
            
            # LSTM Forward
            lstm_out, _ = self.lstm(lstm_input)
            
            # Predict Next Hop using only the LAST time step output
            last_timestep = lstm_out[:, -1, :] # (Batch, Hidden)
            
            out = self.fc(last_timestep)
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
    
    # --- TRANSFER LEARNING ---
    # Load pretrained CNN weights if available
    weights_path = os.path.join('models', 'cnn_weights.pth')
    if os.path.exists(weights_path):
        print(f"Loading pretrained CNN weights from {weights_path}...")
        try:
            # Load state dict
            pretrained_dict = torch.load(weights_path, map_location=device)
            model_dict = model.cnn.state_dict()
            
            # Filter matches just in case
            pretrained_dict = {k: v for k, v in pretrained_dict.items() if k in model_dict}
            
            # Update and load
            model_dict.update(pretrained_dict)
            model.cnn.load_state_dict(model_dict)
            print("✔ Pretrained CNN weights loaded successfully!")
            
            # FREEZE CNN LAYERS
            # This is critical: We want the LSTM to learn to use these perfect features
            # before we try to fine-tune them.
            print("❄ Freezing CNN layers for transfer learning...")
            for param in model.cnn.parameters():
                param.requires_grad = False
                
        except Exception as e:
            print(f"⚠ Failed to load weights: {e}")
    else:
        print("⚠ No pretrained weights found at models/cnn_weights.pth. Training from scratch.")
    
    # --- 6. Train ---
    epochs = 50 # Increased for Colab
    criterion = nn.CrossEntropyLoss()
    
    # AdamW + Cosine Annealing Scheduler
    # Reduced LR for more stable training
    # Optimizer: Only optimize parameters that require gradients
    optimizer = torch.optim.AdamW(filter(lambda p: p.requires_grad, model.parameters()), lr=5e-4, betas=(0.9, 0.999), weight_decay=1e-4)
    scheduler = torch.optim.lr_scheduler.CosineAnnealingLR(optimizer, T_max=epochs, eta_min=1e-6)
    
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
            
            # Gradient clipping for stability
            torch.nn.utils.clip_grad_norm_(model.parameters(), max_norm=1.0)
            
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
        
        # Step scheduler (cosine annealing doesn't need loss)
        scheduler.step()
        
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
