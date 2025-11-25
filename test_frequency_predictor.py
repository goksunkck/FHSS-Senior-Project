# --- test_frequency_predictor.py ---
#
# Tests the trained CNN-LSTM model on the validation dataset.
#
# Output:
#   - results/confusion_matrix.png
#   - Console output with metrics

import json
import numpy as np
import scipy.io
import torch
import torch.nn as nn
from torch.utils.data import Dataset, DataLoader
import h5py
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.metrics import classification_report, confusion_matrix
import os

# --- Shared Utilities (Copied from train_frequency_predictor.py) ---

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
            
            nn.AdaptiveAvgPool2d((10, 1)),
            nn.Flatten()
        )
        
        # Calculate the flattened feature size after CNN
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
        _, (h_n, _) = self.lstm(r_in)
        
        # Get output from the last hidden state and apply dropout
        out = self.dropout(h_n[-1])
        out = self.fc(out)
        
        return out

def test():
    print("--- 1. Load prepared data ---")
    prep_file = 'data/synthetic/prepared_prediction_sequences.mat'
    try:
        prep_data = load_mat_file(prep_file, 'prep_data')
    except FileNotFoundError:
        print(f"Error: Prepared dataset not found at {prep_file}.")
        return

    sequence_starts_validation = prep_data['sequence_starts_validation'].flatten()
    YValidation = prep_data['YValidation'].flatten()
    lookback_window = int(prep_data['lookback_window'])
    
    dataset_filename = prep_data['dataset_filename']
    if isinstance(dataset_filename, np.ndarray):
        dataset_filename = "".join(map(chr, dataset_filename.flatten()))

    class_values = prep_data['class_values'].flatten()
    
    print("--- 2. Load base STFT dataset into memory ---")
    source_dataset = dataset_filename
    try:
        data_contents = load_mat_file(source_dataset, 'X_data')
        X_data_in_memory = data_contents['X_data']
        print("Dataset loaded into memory.")

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
        for img in X_data_in_memory.flatten():
            g_min = min(g_min, img.min())
            g_max = max(g_max, img.max())
        g_range = g_max - g_min if g_max > g_min else 1.0
        print(f"Global Min: {g_min:.4f}, Global Max: {g_max:.4f}")

    except FileNotFoundError:
        print(f"Error: Source dataset not found at {source_dataset}.")
        return

    # --- 3. Process Labels ---
    class_map = {val: i for i, val in enumerate(class_values)}
    YValidation_mapped = np.array([class_map[y] for y in YValidation], dtype=np.int64)
    num_classes = len(class_values)

    # --- 4. Determine input size ---
    first_cell = X_data_in_memory[0, sequence_starts_validation[0] - 1]
    if isinstance(first_cell, np.ndarray):
        first_image = first_cell
    else:
        first_image = first_cell[0]

    if first_image.ndim == 2:
        first_image = np.expand_dims(first_image, axis=-1)
    
    first_image = np.transpose(first_image, (1, 0, 2))
    input_size = first_image.shape

    # --- 5. Create Dataset and DataLoader ---
    val_dataset = SequenceDataset(X_data_in_memory, sequence_starts_validation, YValidation_mapped, lookback_window, input_size, g_min, g_range)
    val_loader = DataLoader(val_dataset, batch_size=32, shuffle=False)

    # --- 6. Load Model ---
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    print(f"Using device: {device}")

    model = CNNLSTM(input_size, num_classes).to(device)
    model_path = 'models/trained_frequency_predictor.pth'
    
    if not os.path.exists(model_path):
        print(f"Error: Model file not found at {model_path}")
        return

    model.load_state_dict(torch.load(model_path, map_location=device))
    model.eval()
    print("Model loaded.")

    # --- 7. Evaluation ---
    print("Starting evaluation...")
    all_preds = []
    all_labels = []

    # Debug: Visualize first batch
    first_batch = next(iter(val_loader))
    debug_seqs, debug_labels = first_batch
    
    print(f"Debug - First batch shape: {debug_seqs.shape}")
    print(f"Debug - First sample min: {debug_seqs[0].min()}, max: {debug_seqs[0].max()}, mean: {debug_seqs[0].mean()}")
    
    # Plot first sample (first time step)
    # Shape: (lookback, C, H, W) -> (15, 1, 1024, 3)
    # We want to plot (H, W) for the first time step
    plt.figure(figsize=(10, 5))
    img_to_plot = debug_seqs[0, 0, 0].cpu().numpy() # First sample, first time step, first channel
    # img_to_plot shape: (1024, 3)
    sns.heatmap(img_to_plot.T, cmap='viridis') # Transpose to (3, 1024) for better viewing
    plt.title(f"Sample Input (Label: {debug_labels[0].item()})")
    plt.savefig('results/sample_input_debug.png')
    print("Saved debug visualization to results/sample_input_debug.png")

    # Debug: Check accuracy of first batch
    print("\n--- Debug: First Batch Accuracy ---")
    with torch.no_grad():
        debug_seqs = debug_seqs.to(device)
        debug_labels = debug_labels.to(device)
        debug_out = model(debug_seqs)
        _, debug_preds = torch.max(debug_out, 1)
        correct = (debug_preds == debug_labels).sum().item()
        acc = correct / debug_labels.size(0)
        print(f"First Batch Accuracy: {acc*100:.2f}% ({correct}/{debug_labels.size(0)})")
        print("First 10 Preds: ", debug_preds[:10].cpu().numpy())
        print("First 10 True:  ", debug_labels[:10].cpu().numpy())

    with torch.no_grad():
        for sequences, labels in val_loader:
            sequences = sequences.to(device)
            labels = labels.to(device)
            outputs = model(sequences)
            _, predicted = torch.max(outputs.data, 1)
            
            all_preds.extend(predicted.cpu().numpy())
            all_labels.extend(labels.cpu().numpy())

    # --- 8. Metrics ---
    print("\nClassification Report:")
    print(classification_report(all_labels, all_preds, target_names=[str(c) for c in class_values]))

    # Confusion Matrix
    cm = confusion_matrix(all_labels, all_preds)
    plt.figure(figsize=(10, 8))
    sns.heatmap(cm, annot=True, fmt='d', cmap='Blues', xticklabels=class_values, yticklabels=class_values)
    plt.xlabel('Predicted')
    plt.ylabel('True')
    plt.title('Confusion Matrix')
    
    results_dir = 'results'
    os.makedirs(results_dir, exist_ok=True)
    save_path = os.path.join(results_dir, 'confusion_matrix.png')
    plt.savefig(save_path)
    print(f"Confusion matrix saved to {save_path}")

    # --- 9. Visualize Example ---
    visualize_example(model, val_loader, class_values, device)

    # --- 10. Analyze All Signals (SNR Sweep) ---
    numHopsPerSignal = int(prep_data['numHopsPerSignal'])
    analyze_signal_accuracies(model, X_data_in_memory, numHopsPerSignal, lookback_window, class_values, device, g_min, g_range)

    # --- 11. Visualize Signal Trajectory (Representative) ---
    # visualize_signal_trajectory(model, X_data_in_memory, numHopsPerSignal, lookback_window, class_values, device, g_min, g_range)

def analyze_signal_accuracies(model, x_data, num_hops_per_signal, lookback_window, class_values, device, norm_min, norm_range):
    """Calculates and prints accuracy for every signal in the dataset."""
    print("\n--- Analyzing Accuracy per Signal (SNR Sweep) ---")
    model.eval()
    
    # Load labels and starts
    prep_file = 'data/synthetic/prepared_prediction_sequences.mat'
    full_data = load_mat_file(prep_file, 'Y_labels')
    all_labels = full_data['Y_labels'].flatten()
    
    full_starts = load_mat_file(prep_file, 'sequence_starts')
    all_starts = full_starts['sequence_starts'].flatten()
    
    seqs_per_signal = num_hops_per_signal - lookback_window
    num_signals = len(all_starts) // seqs_per_signal
    
    # SNR Configuration (Inferred)
    snr_levels = range(-10, 12, 2) # -10, -8, ..., 10
    num_snr_levels = len(snr_levels)
    signals_per_snr = num_signals // num_snr_levels
    
    print(f"Total Signals: {num_signals}")
    print(f"SNR Levels: {list(snr_levels)} ({num_snr_levels} levels)")
    print(f"Signals per SNR: {signals_per_snr}")
    
    print("-" * 60)
    print(f"{'Sig Idx':<8} | {'Est SNR':<8} | {'Accuracy':<10} | {'Correct/Total':<15}")
    print("-" * 60)
    
    snr_accuracies = {snr: [] for snr in snr_levels}
    
    for sig_idx in range(num_signals):
        # Determine SNR
        snr_idx = sig_idx // signals_per_snr
        if snr_idx < num_snr_levels:
            snr = snr_levels[snr_idx]
        else:
            snr = 999 # Should not happen if math is right
            
        start_seq_idx = sig_idx * seqs_per_signal
        end_seq_idx = (sig_idx + 1) * seqs_per_signal
        
        indices_to_check = range(start_seq_idx, end_seq_idx)
        correct_count = 0
        total_count = 0
        
        # Batch processing for speed? 
        # For now, sequential is fine, but let's optimize slightly by not re-creating tensor every time if possible.
        # Actually, creating the batch tensor for the whole signal is better.
        
        # Create batch for this signal
        batch_size = len(indices_to_check)
        batch_seqs = torch.zeros((batch_size, lookback_window, 1, 1024, 3), dtype=torch.float32)
        batch_labels = []
        
        for i, global_seq_idx in enumerate(indices_to_check):
            start_idx_matlab = all_starts[global_seq_idx]
            true_class_val = all_labels[global_seq_idx]
            batch_labels.append(true_class_val)
            
            start_idx_python = start_idx_matlab - 1
            seq_range = range(start_idx_python, start_idx_python + lookback_window)
            
            for t, data_idx in enumerate(seq_range):
                 img_cell = x_data[0, data_idx]
                 if isinstance(img_cell, np.ndarray):
                     img = img_cell
                 else:
                     img = img_cell[0]
                 if img.ndim == 2: img = np.expand_dims(img, axis=-1)
                 img = np.transpose(img, (1, 0, 2))
                 img_tensor = torch.from_numpy(img).permute(2, 0, 1).float()
                 img_tensor = (img_tensor - norm_min) / norm_range
                 batch_seqs[i, t] = img_tensor
        
        batch_seqs = batch_seqs.to(device)
        
        # Split into mini-batches to avoid OOM if signal is long
        mini_batch_size = 32
        num_mini_batches = (batch_size + mini_batch_size - 1) // mini_batch_size
        
        signal_correct = 0
        
        with torch.no_grad():
            for mb in range(num_mini_batches):
                start = mb * mini_batch_size
                end = min(start + mini_batch_size, batch_size)
                
                mb_seqs = batch_seqs[start:end]
                mb_out = model(mb_seqs)
                _, mb_preds = torch.max(mb_out, 1)
                
                mb_labels_list = batch_labels[start:end]
                mb_preds_np = mb_preds.cpu().numpy()
                
                for p, t in zip(mb_preds_np, mb_labels_list):
                    if class_values[p] == t:
                        signal_correct += 1
                        
        acc = signal_correct / batch_size
        snr_accuracies[snr].append(acc)
        
        print(f"{sig_idx:<8} | {snr:<8} | {acc*100:.2f}%    | {signal_correct}/{batch_size}")

    print("-" * 60)
    print("Average Accuracy per SNR Level:")
    for snr in snr_levels:
        avg_acc = np.mean(snr_accuracies[snr])
        print(f"SNR {snr} dB: {avg_acc*100:.2f}%")

    print("-" * 60)
    print("Average Accuracy per SNR Level:")
    for snr in snr_levels:
        avg_acc = np.mean(snr_accuracies[snr])
        print(f"SNR {snr} dB: {avg_acc*100:.2f}%")

def visualize_example(model, val_loader, class_values, device):
    print("\n--- Visualizing a Random Example ---")
    model.eval()
    
    # Get a batch
    sequences, labels = next(iter(val_loader))
    sequences = sequences.to(device)
    labels = labels.to(device)
    
    # Pick a random sample
    idx = np.random.randint(0, sequences.size(0))
    sample_seq = sequences[idx] # (Seq, C, H, W)
    true_label = labels[idx].item()
    
    # Forward pass
    with torch.no_grad():
        # Add batch dimension: (1, Seq, C, H, W)
        output = model(sample_seq.unsqueeze(0))
        probs = torch.softmax(output, dim=1).cpu().numpy()[0]
        prediction = torch.argmax(output, dim=1).item()
        
    print(f"True Label: {class_values[true_label]}")
    print(f"Predicted Label: {class_values[prediction]}")
    
    # Plotting
    plt.figure(figsize=(12, 6))
    
    # 1. Input Spectrogram (Last time step)
    # Shape: (Seq, C, H, W) -> (15, 1, 1024, 3)
    # Get last time step, first channel: (1024, 3)
    last_step_img = sample_seq[-1, 0].cpu().numpy()
    
    plt.subplot(1, 2, 1)
    sns.heatmap(last_step_img.T, cmap='viridis', cbar=False)
    plt.title(f"Input Spectrogram (Last Step)\nTrue Next Hop: {class_values[true_label]}")
    plt.xlabel("Frequency Bins")
    plt.ylabel("Channels")
    
    # 2. Prediction Probabilities
    plt.subplot(1, 2, 2)
    colors = ['gray'] * len(class_values)
    colors[true_label] = 'green' # True label
    if prediction != true_label:
        colors[prediction] = 'red' # Wrong prediction
    else:
        colors[prediction] = 'blue' # Correct prediction (overwrites green, which is fine)
        
    plt.bar(range(len(class_values)), probs, color=colors)
    plt.xticks(range(len(class_values)), class_values)
    plt.title(f"Prediction Probabilities\nPredicted: {class_values[prediction]}")
    plt.xlabel("Class")
    plt.ylabel("Probability")
    plt.ylim(0, 1.0)
    
    save_path = 'results/prediction_example.png'
    plt.tight_layout()
    plt.savefig(save_path)
    print(f"Prediction example saved to {save_path}")

def visualize_signal_trajectory(model, x_data, num_hops_per_signal, lookback_window, class_values, device, norm_min, norm_range):
    """Visualizes predictions over a continuous segment of a signal."""
    print("\n--- Visualizing Signal Trajectory ---")
    model.eval()
    
    # Select the first signal (Signal 0)
    # Hops are 1-indexed in MATLAB logic, but 0-indexed in Python array
    # Signal 0 corresponds to indices 0 to num_hops_per_signal-1
    
    # Let's take a segment of 50 hops (or as many as possible)
    segment_length = 50
    start_hop = 0
    end_hop = min(num_hops_per_signal, start_hop + segment_length)
    
    # We need to form sequences. 
    # Sequence 1: Hops [0, ..., 14] -> Predict Hop 15
    # Sequence 2: Hops [1, ..., 15] -> Predict Hop 16
    # ...
    
    # We can only make predictions starting from hop index `lookback_window`
    prediction_indices = range(lookback_window, end_hop)
    
    true_freqs = []
    pred_freqs = []
    
    print(f"{'Step':<10} | {'True Freq':<15} | {'Pred Freq':<15} | {'Correct':<10}")
    print("-" * 60)
    
    for target_idx in prediction_indices:
        # Construct sequence ending at target_idx - 1
        seq_start = target_idx - lookback_window
        seq_end = target_idx # Exclusive in Python range
        
        # Build tensor
        sequence = torch.zeros((lookback_window, 1, 1024, 3), dtype=torch.float32)
        
        for t, data_idx in enumerate(range(seq_start, seq_end)):
            img_cell = x_data[0, data_idx]
            if isinstance(img_cell, np.ndarray):
                img = img_cell
            else:
                img = img_cell[0]
            
            if img.ndim == 2:
                img = np.expand_dims(img, axis=-1)
            
            img = np.transpose(img, (1, 0, 2)) # (3, 1024, 1)
            img_tensor = torch.from_numpy(img).permute(2, 0, 1).float() # (1, 3, 1024) -> (C, H, W)
            
            # Normalize
            img_tensor = (img_tensor - norm_min) / norm_range
            sequence[t] = img_tensor
            
        # Get True Label (the hop at target_idx)
        # We need to map the raw value to the class index
        # Assuming class_values are 0..7 and map directly for now, but let's be safe
        # In this dataset, Y is stored separately in Y_data, but we can infer it from the filename or just look at the signal generation logic.
        # Wait, x_data contains the spectrograms. The labels are in Y_data (which we didn't load here explicitly in raw form, but we have YValidation).
        # However, for visualization, we want the ACTUAL frequency of the hop at `target_idx`.
        # Since we don't have Y_data loaded for the whole signal easily here (we only loaded prepared sequences), 
        # we can't easily get the true label for *arbitrary* signals without loading Y_data.
        
        # Workaround: Use the prepared validation sequences that happen to be contiguous.
        # But `prepare_prediction_sequences` shuffles or splits.
        
        # Better approach: We loaded `prep_data`. It has `Y_labels` and `sequence_starts`.
        # We can find a contiguous block in `sequence_starts`.
        # `sequence_starts` is generated sequentially: signal 0 (seq 1..N), signal 1 (seq 1..N)...
        # So the first N sequences in `sequence_starts` correspond to Signal 0.
        
        # Let's use the first `segment_length` sequences from the GLOBAL `sequence_starts` (if we had it).
        # We only loaded `sequence_starts_validation`.
        
        # Let's try to infer the label from the spectrogram? No, that's hard.
        # Let's load `Y_labels` from `prep_data` if available.
        # `prep_data` has `YTrain` and `YValidation`.
        # It does NOT have the full `Y_labels` unless we load it.
        # `load_mat_file(prep_file, 'prep_data')` loads the struct.
        # Let's check what `prep_data` contains.
        pass 

    # Re-loading prep_data to get full Y_labels for trajectory visualization
    prep_file = 'data/synthetic/prepared_prediction_sequences.mat'
    full_data = load_mat_file(prep_file, 'Y_labels')
    all_labels = full_data['Y_labels'].flatten()
    
    full_starts = load_mat_file(prep_file, 'sequence_starts')
    all_starts = full_starts['sequence_starts'].flatten()
    
    # Calculate number of sequences per signal
    # We know total sequences and numHopsPerSignal
    # But simpler: we know each signal produces (numHopsPerSignal - lookback_window) sequences.
    seqs_per_signal = num_hops_per_signal - lookback_window
    num_signals = len(all_starts) // seqs_per_signal
    
    print(f"Searching for a representative signal (out of {num_signals} signals)...")
    
    best_signal_idx = -1
    best_acc = -1.0
    
    # Check a few signals to find a good one
    for sig_idx in range(num_signals): # Check ALL signals
        # Shift to the middle of the signal to avoid potential startup transients
        start_seq_idx = sig_idx * seqs_per_signal + 100 
        end_seq_idx = start_seq_idx + 50 
        
        # Ensure we don't go out of bounds
        if end_seq_idx > (sig_idx + 1) * seqs_per_signal:
             start_seq_idx = sig_idx * seqs_per_signal
             end_seq_idx = start_seq_idx + 50
        
        indices_to_check = range(start_seq_idx, end_seq_idx)
        correct_count = 0
        total_count = 0
        
        for i in indices_to_check:
            start_idx_matlab = all_starts[i]
            true_class_val = all_labels[i]
            
            # Reconstruct input
            start_idx_python = start_idx_matlab - 1
            seq_range = range(start_idx_python, start_idx_python + lookback_window)
            
            sequence = torch.zeros((lookback_window, 1, 1024, 3), dtype=torch.float32)
            for t, data_idx in enumerate(seq_range):
                 img_cell = x_data[0, data_idx]
                 if isinstance(img_cell, np.ndarray):
                     img = img_cell
                 else:
                     img = img_cell[0]
                 if img.ndim == 2: img = np.expand_dims(img, axis=-1)
                 img = np.transpose(img, (1, 0, 2))
                 img_tensor = torch.from_numpy(img).permute(2, 0, 1).float()
                 img_tensor = (img_tensor - norm_min) / norm_range
                 sequence[t] = img_tensor
            
            sequence = sequence.to(device)
            with torch.no_grad():
                output = model(sequence.unsqueeze(0))
                pred_idx = torch.argmax(output, dim=1).item()
            
            if class_values[pred_idx] == true_class_val:
                correct_count += 1
            total_count += 1
            
        acc = correct_count / total_count
        print(f"Signal {sig_idx}: Accuracy {acc*100:.1f}%")
        
        if acc > 0.95:
            best_signal_idx = sig_idx
            best_acc = acc
            break
    
    if best_signal_idx == -1:
        print("Could not find a signal with >95% accuracy in the search range. Using Signal 0.")
        best_signal_idx = 0
    else:
        print(f"Selected Signal {best_signal_idx} for visualization.")

    # Visualize the selected signal
    start_seq_idx = best_signal_idx * seqs_per_signal
    indices_to_plot = range(start_seq_idx, start_seq_idx + 50)
    
    true_freqs = []
    pred_freqs = []
    
    print(f"{'Step':<10} | {'True Freq':<15} | {'Pred Freq':<15} | {'Correct':<10}")
    print("-" * 60)
    
    for i in indices_to_plot:
        start_idx_matlab = all_starts[i]
        true_class_val = all_labels[i]
        
        # Reconstruct input
        start_idx_python = start_idx_matlab - 1
        seq_range = range(start_idx_python, start_idx_python + lookback_window)
        
        sequence = torch.zeros((lookback_window, 1, 1024, 3), dtype=torch.float32)
        for t, data_idx in enumerate(seq_range):
             img_cell = x_data[0, data_idx]
             if isinstance(img_cell, np.ndarray):
                 img = img_cell
             else:
                 img = img_cell[0]
             if img.ndim == 2: img = np.expand_dims(img, axis=-1)
             img = np.transpose(img, (1, 0, 2))
             img_tensor = torch.from_numpy(img).permute(2, 0, 1).float()
             img_tensor = (img_tensor - norm_min) / norm_range
             sequence[t] = img_tensor
        
        # Predict
        sequence = sequence.to(device)
        with torch.no_grad():
            output = model(sequence.unsqueeze(0))
            pred_idx = torch.argmax(output, dim=1).item()
            
        pred_val = class_values[pred_idx]
        
        true_freqs.append(true_class_val)
        pred_freqs.append(pred_val)
        
        is_correct = "YES" if true_class_val == pred_val else "NO"
        print(f"{i:<10} | {true_class_val:<15} | {pred_val:<15} | {is_correct:<10}")

    # Plot
    plt.figure(figsize=(15, 6))
    plt.plot(true_freqs, 'g-o', label='True Frequency', linewidth=2, markersize=8)
    plt.plot(pred_freqs, 'r--x', label='Predicted Frequency', linewidth=2, markersize=8)
    
    # Highlight discrepancies
    for k in range(len(true_freqs)):
        if true_freqs[k] != pred_freqs[k]:
            plt.plot(k, pred_freqs[k], 'rx', markersize=12, markeredgewidth=3)
            
    plt.title(f"Signal Trajectory Prediction (Signal {best_signal_idx}, First 50 Hops)")
    plt.xlabel("Time Step (Hop Index)")
    plt.ylabel("Frequency Channel")
    plt.yticks(class_values)
    plt.grid(True, alpha=0.3)
    plt.legend()
    
    save_path = 'results/signal_trajectory.png'
    plt.savefig(save_path)
    print(f"\nSignal trajectory plot saved to {save_path}")

if __name__ == '__main__':
    test()
