
import torch
import torch.nn as nn
import torch.optim as optim
from torch.utils.data import Dataset, DataLoader
import numpy as np
import h5py
import os
import math
import sys

# --- Configuration ---
BATCH_SIZE = 128
EPOCHS = 100
LEARNING_RATE = 1e-3
LOOKBACK_WINDOW = 12
DEVICE = torch.device("cuda" if torch.cuda.is_available() else "cpu")

# --- 1. Dataset (RAM Loading) ---
class RAMDataset(Dataset):
    def __init__(self, h5_filepath, sequence_starts, labels, lookback_window):
        print(f"preload: Loading dataset from {h5_filepath} into RAM...")
        self.sequence_starts = sequence_starts
        self.labels = labels
        self.lookback_window = lookback_window
        
        with h5py.File(h5_filepath, 'r') as f:
            x_data = f['X_data'][:] 
            
            # Normalize Per Sample
            print("Normalizing...")
            flat = x_data.reshape(x_data.shape[0], -1)
            mins = flat.min(axis=1, keepdims=True)
            maxs = flat.max(axis=1, keepdims=True)
            rngs = maxs - mins
            rngs[rngs == 0] = 1.0
            
            mins = mins.reshape(-1, 1, 1)
            rngs = rngs.reshape(-1, 1, 1)
            
            x_data = (x_data - mins) / rngs
            
            self.data = torch.from_numpy(x_data).float().unsqueeze(1) # (N, 1, 1024, 3)
            print(f"Loaded {self.data.shape}")

    def __len__(self):
        return len(self.sequence_starts)

    def __getitem__(self, idx):
        start = int(self.sequence_starts[idx])
        end = start + self.lookback_window
        
        seq = self.data[start:end] 
        label = self.labels[idx]
        
        return seq, torch.tensor(label, dtype=torch.long)

# ... (Previous Model Code Omitted) ...

def main():
    print("--- Training CNN-Transformer on M-SEQUENCES ---")
    
    # Paths M-SEQUENCE SPECIFIC
    raw_h5 = os.path.join('data', 'synthetic', 'm_sequence_dataset.h5')
    prep_h5 = os.path.join('data', 'synthetic', 'prepared_m_sequence_indices.h5')
    
    if not os.path.exists(raw_h5) or not os.path.exists(prep_h5):
        print(f"Data files missing!\nRaw: {raw_h5}\nPrep: {prep_h5}")
        return
        
    # Read Meta
    with h5py.File(prep_h5, 'r') as f:
        starts_train = f['sequence_starts_train'][:]
        y_train = f['YTrain'][:]
        starts_val = f['sequence_starts_validation'][:]
        y_val = f['YValidation'][:]
        
        # Check if saved meta matches config
        if 'lookback_window' in f.attrs:
            saved_window = f.attrs['lookback_window']
            if saved_window != LOOKBACK_WINDOW:
                print(f"WARNING: Config LOOKBACK_WINDOW={LOOKBACK_WINDOW} but prepared data has {saved_window}!")
        
    print(f"Meta Loaded. Train: {len(starts_train)}, Val: {len(starts_val)}")
    
    # Init Datasets
    train_ds = RAMDataset(raw_h5, starts_train, y_train, LOOKBACK_WINDOW)
    val_ds = RAMDataset(raw_h5, starts_val, y_val, LOOKBACK_WINDOW)
    
    train_loader = DataLoader(train_ds, batch_size=BATCH_SIZE, shuffle=True)
    val_loader = DataLoader(val_ds, batch_size=BATCH_SIZE, shuffle=False)
    
    # Setup Model
    model = CNNTransformer(8).to(DEVICE)
    
    # Load Weights (Transfer Learning from Gold Code CNN is valid because physics is same)
    w_path = os.path.join('models', 'cnn_weights.pth')
    if os.path.exists(w_path):
        print("Loading Pretrained CNN...")
        d = torch.load(w_path, map_location=DEVICE)
        model.cnn.load_state_dict(d, strict=False)
    
    optimizer = optim.AdamW(model.parameters(), lr=LEARNING_RATE)
    scheduler = optim.lr_scheduler.CosineAnnealingLR(optimizer, T_max=EPOCHS)
    criterion = nn.CrossEntropyLoss()
    
    # History Tracking
    history = {'train_loss': [], 'val_loss': [], 'train_acc': [], 'val_acc': []}
    best_val_loss = float('inf')
    
    for epoch in range(EPOCHS):
        model.train()
        total_loss = 0
        correct = 0
        total = 0
        
        for X, y in train_loader:
            X, y = X.to(DEVICE), y.to(DEVICE)
            y = y.view(-1)
            
            optimizer.zero_grad()
            out = model(X)
            loss = criterion(out, y)
            loss.backward()
            optimizer.step()
            
            total_loss += loss.item()
            pred = out.argmax(dim=1)
            correct += (pred == y).sum().item()
            total += y.size(0)
            
        # Validation
        val_loss = 0
        model.eval()
        v_corr = 0
        v_tot = 0
        with torch.no_grad():
            for X, y in val_loader:
                X, y = X.to(DEVICE), y.to(DEVICE)
                y = y.view(-1)
                
                out = model(X)
                v_loss_batch = criterion(out, y)
                val_loss += v_loss_batch.item()
                
                v_corr += (out.argmax(dim=1) == y).sum().item()
                v_tot += y.size(0)
        
        # Metrics
        avg_train_loss = total_loss / len(train_loader)
        avg_val_loss = val_loss / len(val_loader)
        train_acc = 100 * correct / total
        val_acc = 100 * v_corr / v_tot
        
        history['train_loss'].append(avg_train_loss)
        history['val_loss'].append(avg_val_loss)
        history['train_acc'].append(train_acc)
        history['val_acc'].append(val_acc)
        
        scheduler.step()
        
        print(f"Epoch {epoch+1}/{EPOCHS} | Train Loss: {avg_train_loss:.4f} | Val Loss: {avg_val_loss:.4f} | Train Acc: {train_acc:.1f}% | Val Acc: {val_acc:.1f}%")
        
        # Save Best Model (Loss Based)
        if avg_val_loss < best_val_loss:
            best_val_loss = avg_val_loss
            torch.save(model.state_dict(), os.path.join('models', 'm_sequence_transformer_best.pth'))
            print("  --> Best Model Saved (Improved Val Loss)")
            
    # Plotting
    import matplotlib.pyplot as plt
    
    plt.figure(figsize=(12, 5))
    
    plt.subplot(1, 2, 1)
    plt.plot(history['train_loss'], label='Train Loss')
    plt.plot(history['val_loss'], label='Val Loss')
    plt.title('Loss History')
    plt.xlabel('Epoch')
    plt.legend()
    
    plt.subplot(1, 2, 2)
    plt.plot(history['train_acc'], label='Train Acc')
    plt.plot(history['val_acc'], label='Val Acc')
    plt.title('Accuracy History')
    plt.xlabel('Epoch')
    plt.legend()
    
    plt.tight_layout()
    plt.savefig('results/training_history_mseq.png')
    print("Training Complete. History plot saved to results/training_history_mseq.png")

if __name__ == "__main__":
    main()
