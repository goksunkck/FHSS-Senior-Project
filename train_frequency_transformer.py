
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
BATCH_SIZE = 64
EPOCHS = 30
LEARNING_RATE = 1e-3
DEVICE = torch.device("cuda" if torch.cuda.is_available() else "cpu")

# --- 1. Dataset (RAM Loading - Same as Colab Script) ---
class RAMDataset(Dataset):
    def __init__(self, h5_filepath, sequence_starts, labels, lookback_window):
        print(f"preload: Loading dataset from {h5_filepath} into RAM...")
        self.sequence_starts = sequence_starts
        self.labels = labels
        self.lookback_window = lookback_window
        
        with h5py.File(h5_filepath, 'r') as f:
            x_data = f['X_data'][:] # Load ALL into RAM (Fast)
            
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
        
        seq = self.data[start:end] # Slice tensor
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
        pe = pe.unsqueeze(0).transpose(0, 1) # (S, 1, E)
        self.register_buffer('pe', pe)

    def forward(self, x):
        return x + self.pe[:x.size(0), :]

# --- 3. Model Architecture (CNN-Transformer) ---
class CNNTransformer(nn.Module):
    def __init__(self, num_classes, d_model=64, nhead=4, num_layers=2):
        super(CNNTransformer, self).__init__()
        
        # Vision (Must match pre-trained)
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
        self.cnn_fc = nn.Linear(32 * 64, num_classes)
        
        # Symbolic
        self.embedding = nn.Embedding(num_classes, d_model)
        self.pos_encoder = PositionalEncoding(d_model)
        
        # Transformer
        encoder_layers = nn.TransformerEncoderLayer(d_model=d_model, nhead=nhead, dim_feedforward=128, dropout=0.0)
        self.transformer_encoder = nn.TransformerEncoder(encoder_layers, num_layers=num_layers)
        
        self.fc = nn.Linear(d_model, num_classes)

    def forward(self, x):
        # x: (B, S, 1, F, T)
        b, s, c, h, w = x.size()
        c_in = x.view(b * s, c, h, w)
        
        with torch.no_grad():
            c_out = self.cnn(c_in)
            c_flat = c_out.view(c_out.size(0), -1)
            logits = self.cnn_fc(c_flat)
            _, ids = torch.max(logits, dim=1)
            
        emb = self.embedding(ids)
        src = emb.view(b, s, -1).permute(1, 0, 2) # (S, B, E)
        src = self.pos_encoder(src)
        out = self.transformer_encoder(src)
        last = out[-1, :, :]
        return self.fc(last)

def main():
    print("--- Training CNN-Transformer (RAM Mode) ---")
    
    # Paths
    raw_h5 = os.path.join('data', 'synthetic', 'classification_dataset_stft_random_deg10_python.h5')
    prep_h5 = os.path.join('data', 'synthetic', 'prepared_prediction_sequences_python.h5')
    
    if not os.path.exists(raw_h5) or not os.path.exists(prep_h5):
        print("Data files missing!")
        return
        
    # Read Meta
    with h5py.File(prep_h5, 'r') as f:
        starts_train = f['sequence_starts_train'][:]
        y_train = f['YTrain'][:]
        starts_val = f['sequence_starts_validation'][:]
        y_val = f['YValidation'][:]
        
    print(f"Meta Loaded. Train: {len(starts_train)}, Val: {len(starts_val)}")
    
    # Init Datasets
    train_ds = RAMDataset(raw_h5, starts_train, y_train, 12)
    val_ds = RAMDataset(raw_h5, starts_val, y_val, 12)
    
    train_loader = DataLoader(train_ds, batch_size=BATCH_SIZE, shuffle=True)
    val_loader = DataLoader(val_ds, batch_size=BATCH_SIZE, shuffle=False)
    
    # Setup Model
    model = CNNTransformer(8).to(DEVICE)
    
    # Load Weights
    w_path = os.path.join('models', 'cnn_weights.pth')
    if os.path.exists(w_path):
        print("Loading Pretrained CNN...")
        d = torch.load(w_path, map_location=DEVICE)
        model.cnn.load_state_dict(d, strict=False)
    
    optimizer = optim.AdamW(model.parameters(), lr=LEARNING_RATE)
    criterion = nn.CrossEntropyLoss()
    
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
        
        avg_val_loss = val_loss / len(val_loader)
        print(f"Epoch {epoch+1} | Train Loss: {total_loss/len(train_loader):.4f} | Val Loss: {avg_val_loss:.4f} | Train Acc: {100*correct/total:.1f}% | Val Acc: {100*v_corr/v_tot:.1f}%")
        
        torch.save(model.state_dict(), os.path.join('models', 'cnn_transformer.pth'))

if __name__ == "__main__":
    main()
