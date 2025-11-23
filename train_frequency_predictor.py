# --- train_frequency_predictor.py ---
#
# Trains a CNN-LSTM to predict the next hop frequency index using PyTorch.
# This script is a Python conversion of the original MATLAB script.
#
# Requires:
#   - data/synthetic/prepared_prediction_sequences.mat (index + label sets)
#   - data/synthetic/classification_dataset_stft_awgn.mat (STFT tensors)
#
# Output:
#   - models/trained_frequency_predictor.pth (PyTorch model state dictionary)
#   - models/model_metadata.json (training metadata)

import json
import numpy as np
import scipy.io
import torch
import torch.nn as nn
from torch.utils.data import Dataset, DataLoader
import h5py


def load_mat_file(filepath, variable_name):
    """Loads a specific variable from a .mat file."""
    try:
        return scipy.io.loadmat(filepath)
    except NotImplementedError:
        # If scipy.io.loadmat fails, it might be a v7.3 .mat file.
        # h5py is needed for this.
        import h5py
        with h5py.File(filepath, 'r') as f:
            data = {}
            for k, v in f.items():
                if isinstance(v, h5py.Dataset): # Check if it's a dataset
                    # Dereference and copy data to make it mutable
                    value = v[()]
                    # Handle MATLAB's tendency to save everything as arrays
                    if isinstance(value, np.ndarray):
                        if value.shape == (1, 1):
                            data[k] = value[0, 0] # Unpack scalar
                        else:
                            data[k] = value
                    else:
                        data[k] = value
            return data

def train():
    """Main training function."""
    print("--- 1. Load prepared data ---")
    prep_file = 'data/synthetic/prepared_prediction_sequences.mat'
    try:
        prep_data = load_mat_file(prep_file, 'prep_data')
    except FileNotFoundError:
        print(f"Error: Prepared dataset not found at {prep_file}.")
        print("Please run the MATLAB script 'prepare_prediction_sequences.m' first.")
        return

    sequence_starts_train = prep_data['sequence_starts_train'].flatten()
    sequence_starts_validation = prep_data['sequence_starts_validation'].flatten()
    YTrain = prep_data['YTrain'].flatten()
    YValidation = prep_data['YValidation'].flatten()
    lookback_window = int(prep_data['lookback_window'])
    
    dataset_filename = prep_data['dataset_filename']
    if isinstance(dataset_filename, np.ndarray):
        dataset_filename = "".join(map(chr, dataset_filename.flatten()))

    class_values = prep_data['class_values'].flatten()
    
    class_names = prep_data['class_names']


    print("--- 2. Load base STFT dataset into memory ---")
    source_dataset = dataset_filename
    try:
        # Note: loading v7.3 mat files requires h5py.
        data_contents = load_mat_file(source_dataset, 'X_data')
        X_data_in_memory = data_contents['X_data']
        print("Dataset loaded into memory.")

        # If dataset is v7.3, X_data might be an array of object references. We must resolve them.
        if X_data_in_memory.dtype == h5py.ref_dtype:
            print("Resolving HDF5 object references for X_data...")
            resolved_x_data = np.empty(X_data_in_memory.shape, dtype=object)
            with h5py.File(source_dataset, 'r') as f:
                for i in range(X_data_in_memory.size):
                    ref = X_data_in_memory.flat[i]
                    resolved_x_data.flat[i] = f[ref][()]
            X_data_in_memory = resolved_x_data
            print("References resolved.")

        print("Computing normalization stats over the dataset...")
        g_min = float('inf')
        g_max = float('-inf')
        # This loop is memory-efficient for large datasets
        for img in X_data_in_memory.flatten():
            g_min = min(g_min, img.min())
            g_max = max(g_max, img.max())
        g_range = g_max - g_min if g_max > g_min else 1.0
        print(f"Global Min: {g_min:.4f}, Global Max: {g_max:.4f}")

    except FileNotFoundError:
        print(f"Error: Source dataset not found at {source_dataset}.")
        print("Please re-run the MATLAB script 'generate_prediction_dataset.m'.")
        return
    
    print(f"Loaded {len(sequence_starts_train)} training and {len(sequence_starts_validation)} validation sequence indices.")

    # --- 3. Process Labels ---
    # In PyTorch, labels are typically 0-indexed integers for CrossEntropyLoss.
    # We will map the MATLAB class_values to 0, 1, 2, ...
    # And we assume YTrain/YValidation contain the actual class values, not indices.
    
    # Create a mapping from the original class value to a 0-indexed label
    class_map = {val: i for i, val in enumerate(class_values)}
    
    YTrain_mapped = np.array([class_map[y] for y in YTrain], dtype=np.int64)
    YValidation_mapped = np.array([class_map[y] for y in YValidation], dtype=np.int64)

    num_classes = len(class_values)

    print("\n--- Label Inspection ---")
    print(f"Class values (from .mat file): {class_values}")
    print(f"Number of distinct class values: {len(class_values)}")
    print(f"Num classes (for model): {num_classes}")
    print(f"YTrain_mapped type: {YTrain_mapped.dtype}, shape: {YTrain_mapped.shape}")
    print(f"YTrain_mapped unique values: {np.unique(YTrain_mapped)}")
    print(f"First 10 YTrain_mapped labels: {YTrain_mapped[:10]}")
    print("------------------------\n")

    # --- 4. Determine input size ---
    # We need to inspect the first item to get the dimensions
    first_cell = X_data_in_memory[0, sequence_starts_train[0] - 1] # MATLAB is 1-indexed, and data is row-major
    # When loading from a v7.3 MAT file, it might already be a numpy array
    if isinstance(first_cell, np.ndarray):
        first_image = first_cell
    else: # Fallback for standard .mat loading
        first_image = first_cell[0]

    # Ensure the image is 3D (H, W, C) for the model
    if first_image.ndim == 2:
        first_image = np.expand_dims(first_image, axis=-1)
    
    # Transpose from (3, 1024, 1) to (1024, 3, 1) to match model architecture
    first_image = np.transpose(first_image, (1, 0, 2))

    input_size = first_image.shape # (H, W, C)
    print(f'Input STFT size: {input_size}')
    print(f'Num classes: {num_classes}')

    # --- 5. Create PyTorch Datasets and DataLoaders ---
    print("--- 5. Creating Datasets and DataLoaders ---")

    class SequenceDataset(Dataset):
        """Custom Dataset for loading sequences."""
        def __init__(self, x_data, sequence_starts, labels, lookback_window, input_size, norm_min, norm_range):
            self.x_data = x_data
            self.sequence_starts = sequence_starts
            self.labels = labels
            self.lookback_window = lookback_window
            self.input_size = input_size # H, W, C
            self.norm_min = norm_min
            self.norm_range = norm_range

        def __len__(self):
            return len(self.sequence_starts)

        def __getitem__(self, idx):
            start_idx_matlab = self.sequence_starts[idx]
            label = self.labels[idx]
            
            # Convert from MATLAB 1-based index to Python 0-based index
            start_idx_python = start_idx_matlab - 1
            
            seq_range = range(start_idx_python, start_idx_python + self.lookback_window)
            
            # Pre-allocate sequence tensor: (lookback, C, H, W) for PyTorch
            # Input size is (H, W, C), so we permute
            sequence = torch.zeros((self.lookback_window, self.input_size[2], self.input_size[0], self.input_size[1]), dtype=torch.float32)

            for t, data_idx in enumerate(seq_range):
                # Data is stored as a numpy object array, loaded as a row-vector (1, N)
                # Access column `data_idx` from the first row
                img_cell = self.x_data[0, data_idx]
                
                # Ensure we have a numpy array
                if isinstance(img_cell, np.ndarray):
                    img = img_cell
                else:
                    img = img_cell[0]

                # Add channel dimension if it's 2D
                if img.ndim == 2:
                    img = np.expand_dims(img, axis=-1)
                
                # Transpose from (3, 1024, 1) to (1024, 3, 1) to match model architecture
                img = np.transpose(img, (1, 0, 2))
                
                # Transpose from (H, W, C) to (C, H, W) for PyTorch
                img_tensor = torch.from_numpy(img).permute(2, 0, 1).float()
                
                # Normalize the tensor
                img_tensor = (img_tensor - self.norm_min) / self.norm_range
                sequence[t] = img_tensor
            
            return sequence, torch.tensor(label, dtype=torch.long)

    # Instantiate datasets
    train_dataset = SequenceDataset(X_data_in_memory, sequence_starts_train, YTrain_mapped, lookback_window, input_size, g_min, g_range)
    val_dataset = SequenceDataset(X_data_in_memory, sequence_starts_validation, YValidation_mapped, lookback_window, input_size, g_min, g_range)

    # Instantiate dataloaders
    batch_size = 4
    train_loader = DataLoader(train_dataset, batch_size=batch_size, shuffle=True)
    val_loader = DataLoader(val_dataset, batch_size=batch_size, shuffle=False)

    print("--- 6. Define CNN-LSTM architecture ---")

    class CNNLSTM(nn.Module):
        def __init__(self, input_size, num_classes):
            super(CNNLSTM, self).__init__()
            # CNN part - designed to accept (C, H, W)
            self.cnn = nn.Sequential(
                nn.Conv2d(input_size[2], 16, kernel_size=3, padding='same'),
                nn.BatchNorm2d(16),
                nn.ReLU(),
                nn.MaxPool2d(kernel_size=(2, 1), stride=(2, 1)),

                nn.Conv2d(16, 32, kernel_size=3, padding='same'),
                nn.BatchNorm2d(32),
                nn.ReLU(),
                nn.MaxPool2d(kernel_size=(2, 1), stride=(2, 1)),
                
                nn.AdaptiveAvgPool2d((10, 1)), # Adding this layer
                nn.Flatten()
            )
            
            # Calculate the flattened feature size after CNN
            # We pass a dummy tensor through the CNN to find the output size
            with torch.no_grad():
                dummy_input = torch.zeros(1, input_size[2], input_size[0], input_size[1])
                cnn_out_size = self.cnn(dummy_input).shape[1]

            # LSTM part
            self.lstm = nn.LSTM(input_size=cnn_out_size, hidden_size=100, batch_first=True)
            self.dropout = nn.Dropout(0.5)
            
            # Fully connected part
            self.fc = nn.Linear(100, num_classes)

        def forward(self, x):
            # x shape: (batch, lookback, C, H, W)
            batch_size, seq_len, C, H, W = x.size()
            
            # Reshape to pass each time step through CNN
            c_in = x.view(batch_size * seq_len, C, H, W)
            
            # CNN features
            c_out = self.cnn(c_in)
            
            # Reshape back for LSTM: (batch, lookback, features)
            r_in = c_out.view(batch_size, seq_len, -1)
            
            # LSTM
            # We only need the last hidden state
            _, (h_n, _) = self.lstm(r_in)
            
            # Get output from the last hidden state and apply dropout
            out = self.dropout(h_n[-1])
            out = self.fc(out)
            
            return out

    # PyTorch device setup
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    print(f"Using device: {device}")

    # Instantiate the model
    model = CNNLSTM(input_size, num_classes).to(device)
    print(model)

    # --- 7. Training options ---
    epochs = 20
    learning_rate = 1e-4
    criterion = nn.CrossEntropyLoss()
    optimizer = torch.optim.Adam(model.parameters(), lr=learning_rate)
    
    # --- 8. Train ---
    print('--- Starting training... ---')
    for epoch in range(epochs):
        model.train()
        running_loss = 0.0
        for i, (sequences, labels) in enumerate(train_loader):
            sequences = sequences.to(device)
            labels = labels.to(device)
            
            # Zero gradients
            optimizer.zero_grad()
            
            # Forward pass
            outputs = model(sequences)
            loss = criterion(outputs, labels)
            
            # Backward and optimize
            loss.backward()
            optimizer.step()
            
            running_loss += loss.item()
            if (i + 1) % 100 == 0:
                print(f'Epoch [{epoch+1}/{epochs}], Step [{i+1}/{len(train_loader)}], Loss: {loss.item():.4f}')

        # Validation
        model.eval()
        with torch.no_grad():
            correct = 0
            total = 0
            val_loss = 0
            for sequences, labels in val_loader:
                sequences = sequences.to(device)
                labels = labels.to(device)
                outputs = model(sequences)
                val_loss += criterion(outputs, labels).item()
                _, predicted = torch.max(outputs.data, 1)
                total += labels.size(0)
                correct += (predicted == labels).sum().item()
        
        avg_val_loss = val_loss / len(val_loader)
        val_accuracy = 100 * correct / total
        print(f'Epoch [{epoch+1}/{epochs}] - Validation Loss: {avg_val_loss:.4f}, Validation Accuracy: {val_accuracy:.2f}%')

    print('--- Training complete. ---')


if __name__ == '__main__':
    train()
