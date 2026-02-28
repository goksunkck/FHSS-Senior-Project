
import torch
import torch.nn as nn
import torch.optim as optim
from torch.utils.data import Dataset, DataLoader
import numpy as np
import h5py
import os
import math
import sys
import matplotlib.pyplot as plt

# --- Configuration ---
BATCH_SIZE = 128
EPOCHS = 50
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

# --- 2. Positional Encoding ---
class PositionalEncoding(nn.Module):
    def __init__(self, d_model, max_len=5000):
        super(PositionalEncoding, self).__init__()
        pe = torch.zeros(max_len, d_model)
        position = torch.arange(0, max_len, dtype=torch.float).unsqueeze(1)
        div_term = torch.exp(torch.arange(0, d_model, 2).float() * (-math.log(10000.0) / d_model))
        pe[:, 0::2] = torch.sin(position * div_term)
        pe[:, 1::2] = torch.cos(position * div_term)
        pe = pe.unsqueeze(0).transpose(0, 1) 
        self.register_buffer('pe', pe)

    def forward(self, x):
        return x + self.pe[:x.size(0), :]

# --- 3. Model Architecture (CNN-Transformer) ---
class CNNTransformer(nn.Module):
    def __init__(self, num_classes, d_model=64, nhead=4, num_layers=2):
        super(CNNTransformer, self).__init__()
        
        # Vision (Must match train_cnn.py EXACTLY)
        self.cnn = nn.Sequential(
            # Layer 1
            nn.Conv2d(1, 32, kernel_size=(5, 1), stride=1, padding=(2, 0)),
            nn.InstanceNorm2d(32),
            nn.ReLU(),
            nn.MaxPool2d(kernel_size=(4, 1), stride=(4, 1)),
            
            # Layer 2
            nn.Conv2d(32, 64, kernel_size=(5, 1), stride=1, padding=(2, 0)),
            nn.InstanceNorm2d(64),
            nn.ReLU(),
            nn.MaxPool2d(kernel_size=(4, 1), stride=(4, 1)),

            # Layer 3
            nn.Conv2d(64, 32, kernel_size=(5, 1), stride=1, padding=(2, 0)),
            nn.InstanceNorm2d(32),
            nn.ReLU(),
            
            # Layer 4
            nn.Conv2d(32, 32, kernel_size=(5, 1), stride=1, padding=(2, 0)),
            nn.InstanceNorm2d(32),
            nn.ReLU()
        )
        
        self.flat_size = 32 * 64 * 3  # 32 (Channels) * 64 (Height: 1024->256->64) * 3 (Width: Preserved) 
        self.cnn_fc = nn.Linear(self.flat_size, num_classes)
        
        # Symbolic
        self.embedding = nn.Embedding(num_classes, d_model)
        self.pos_encoder = PositionalEncoding(d_model)
        
        # Transformer
        encoder_layers = nn.TransformerEncoderLayer(d_model=d_model, nhead=nhead, dim_feedforward=128, dropout=0.0)
        self.transformer_encoder = nn.TransformerEncoder(encoder_layers, num_layers=num_layers)
        
        # We output 3 bits instead of a single class integer (8 classes)
        self.fc = nn.Linear(d_model, 3)

    def forward(self, x):
        b, s, c, h, w = x.size()
        c_in = x.view(b * s, c, h, w)
        
        with torch.no_grad():
            c_out = self.cnn(c_in)
            c_flat = c_out.view(c_out.size(0), -1)
            logits = self.cnn_fc(c_flat)
            _, ids = torch.max(logits, dim=1)
            
        emb = self.embedding(ids)
        src = emb.view(b, s, -1).permute(1, 0, 2) 
        src = self.pos_encoder(src)
        out = self.transformer_encoder(src)
        last = out[-1, :, :]
        return self.fc(last)

def main():
    print("--- Training CNN-Transformer on GOLD CODES ---")
    
    # Paths (Gold Code Specific)
    raw_h5 = os.path.join('data', 'synthetic', 'classification_dataset_stft_random_deg10_python.h5')
    prep_h5 = os.path.join('data', 'synthetic', 'prepared_prediction_sequences_python.h5')
    
    if not os.path.exists(raw_h5) or not os.path.exists(prep_h5):
        print(f"Data files missing!\nRaw: {raw_h5}\nPrep: {prep_h5}")
        return
        
    # Read Meta
    with h5py.File(prep_h5, 'r') as f:
        starts_train = f['sequence_starts_train'][:]
        y_train = f['YTrain'][:]
        starts_val = f['sequence_starts_validation'][:]
        y_val = f['YValidation'][:]
        
    print(f"Meta Loaded. Train: {len(starts_train)}, Val: {len(starts_val)}")
    
    # Init Datasets
    train_ds = RAMDataset(raw_h5, starts_train, y_train, LOOKBACK_WINDOW)
    val_ds = RAMDataset(raw_h5, starts_val, y_val, LOOKBACK_WINDOW)
    
    train_loader = DataLoader(train_ds, batch_size=BATCH_SIZE, shuffle=True)
    val_loader = DataLoader(val_ds, batch_size=BATCH_SIZE, shuffle=False)
    
    # Setup Model
    model = CNNTransformer(8).to(DEVICE)
    
    # Load Weights (Transfer Learning)
    w_path = os.path.join('models', 'cnn_weights.pth')
    if os.path.exists(w_path):
        print("Loading Pretrained CNN (Full)...")
        d = torch.load(w_path, map_location=DEVICE)
        
        new_cnn_dict = {}
        new_fc_dict = {}
        
        for k, v in d.items():
            if k.startswith('cnn.'):
                name = k[4:] 
                new_cnn_dict[name] = v
            elif k.startswith('fc.'):
                name = k[3:]
                new_fc_dict[name] = v
                
        model.cnn.load_state_dict(new_cnn_dict, strict=True)
        model.cnn_fc.weight.data = new_fc_dict['weight']
        model.cnn_fc.bias.data = new_fc_dict['bias']
        print("Loaded CNN Backbone + Classification Head successfully.")
    else:
        print("WARNING: Pretrained weights not found!")
    
    optimizer = optim.AdamW(model.parameters(), lr=LEARNING_RATE)
    scheduler = optim.lr_scheduler.CosineAnnealingLR(optimizer, T_max=EPOCHS)
    # We use BCEWithLogitsLoss because we are predicting 3 independent bits
    criterion = nn.BCEWithLogitsLoss()
    
    # History Tracking
    history = {'train_loss': [], 'val_loss': [], 'train_acc': [], 'val_acc': [], 'train_seq_acc': [], 'val_seq_acc': []}
    best_val_loss = float('inf')
    
    for epoch in range(EPOCHS):
        model.train()
        total_loss = 0
        correct = 0
        total = 0
        correct_seq = 0
        total_seq = 0
        
        for X, y in train_loader:
            X, y = X.to(DEVICE), y.to(DEVICE)
            y = y.view(-1)
            
            # Convert class label (0-7) to 3-bit binary
            y_bits = torch.zeros(y.size(0), 3).to(DEVICE)
            for i in range(3):
                y_bits[:, i] = (y >> (2 - i)) & 1
            
            optimizer.zero_grad()
            out = model(X)
            # out is [batch, 3] and y_bits is [batch, 3]
            loss = criterion(out, y_bits)
            loss.backward()
            optimizer.step()
            
            total_loss += loss.item()
            
            pred_bits = (out > 0).float()
            correct += (pred_bits == y_bits).sum().item()
            total += y.size(0) * 3  # 3 bits per sample
            
            correct_seq += (pred_bits == y_bits).all(dim=1).sum().item()
            total_seq += y.size(0)
            
        # Validation
        val_loss = 0
        model.eval()
        v_corr = 0
        v_tot = 0
        v_seq_corr = 0
        v_seq_tot = 0
        with torch.no_grad():
            for X, y in val_loader:
                X, y = X.to(DEVICE), y.to(DEVICE)
                y = y.view(-1)
                
                y_bits = torch.zeros(y.size(0), 3).to(DEVICE)
                for i in range(3):
                    y_bits[:, i] = (y >> (2 - i)) & 1
                
                out = model(X)
                v_loss_batch = criterion(out, y_bits)
                val_loss += v_loss_batch.item()
                
                pred_bits = (out > 0).float()
                v_corr += (pred_bits == y_bits).sum().item()
                v_tot += y.size(0) * 3
                
                v_seq_corr += (pred_bits == y_bits).all(dim=1).sum().item()
                v_seq_tot += y.size(0)
        
        # Metrics
        avg_train_loss = total_loss / len(train_loader)
        avg_val_loss = val_loss / len(val_loader)
        train_acc = 100 * correct / total
        val_acc = 100 * v_corr / v_tot
        train_seq_acc = 100 * correct_seq / total_seq
        val_seq_acc = 100 * v_seq_corr / v_seq_tot
        
        history['train_loss'].append(avg_train_loss)
        history['val_loss'].append(avg_val_loss)
        history['train_acc'].append(train_acc)
        history['val_acc'].append(val_acc)
        history['train_seq_acc'].append(train_seq_acc)
        history['val_seq_acc'].append(val_seq_acc)
        
        scheduler.step()
        
        print(f"Epoch {epoch+1}/{EPOCHS} | Train Loss: {avg_train_loss:.4f} | Val Loss: {avg_val_loss:.4f} | Val Bit Acc: {val_acc:.1f}% | Val Hop Acc: {val_seq_acc:.1f}%")
        
        if avg_val_loss < best_val_loss:
            best_val_loss = avg_val_loss
            torch.save(model.state_dict(), os.path.join('models', 'gold_code_transformer_best.pth'))
            print("  --> Best Model Saved")
            
    # Plotting
    plt.figure(figsize=(12, 5))
    
    plt.subplot(1, 2, 1)
    plt.plot(history['train_loss'], label='Train Loss')
    plt.plot(history['val_loss'], label='Val Loss')
    plt.title('Gold Code Loss History')
    plt.xlabel('Epoch')
    plt.legend()
    
    plt.subplot(1, 2, 2)
    plt.plot(history['train_acc'], label='Train Acc')
    plt.plot(history['val_acc'], label='Val Acc')
    plt.title('Gold Code Accuracy History')
    plt.xlabel('Epoch')
    plt.legend()
    
    plt.tight_layout()
    plt.savefig('results/training_history_gold.png')
    print("Training Complete. History plot saved to results/training_history_gold.png")

if __name__ == "__main__":
    main()
