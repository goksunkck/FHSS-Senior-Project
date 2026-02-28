
import torch
import torch.nn as nn
from torch.utils.data import Dataset, DataLoader
import numpy as np
import h5py
import os
import math
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.metrics import confusion_matrix

# --- Config ---
BATCH_SIZE = 128
LOOKBACK_WINDOW = 12
DEVICE = torch.device("cuda" if torch.cuda.is_available() else "cpu")
MODEL_PATH = os.path.join('models', 'gold_code_transformer_best.pth')

# --- 1. Re-Define Architecture (Matches train_frequency_transformer.py) ---
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

class CNNTransformer(nn.Module):
    def __init__(self, num_classes, d_model=64, nhead=4, num_layers=2):
        super(CNNTransformer, self).__init__()
        
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
        
        self.flat_size = 32 * 64 * 3 
        self.cnn_fc = nn.Linear(self.flat_size, num_classes)
        
        self.embedding = nn.Embedding(num_classes, d_model)
        self.pos_encoder = PositionalEncoding(d_model)
        
        encoder_layers = nn.TransformerEncoderLayer(d_model=d_model, nhead=nhead, dim_feedforward=128, dropout=0.0)
        self.transformer_encoder = nn.TransformerEncoder(encoder_layers, num_layers=num_layers)
        
        self.fc = nn.Linear(d_model, num_classes)

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

# --- 2. Dataset ---
class RAMDataset(Dataset):
    def __init__(self, h5_filepath, sequence_starts, labels, lookback_window):
        print(f"Loading dataset from {h5_filepath}...")
        self.sequence_starts = sequence_starts
        self.labels = labels
        self.lookback_window = lookback_window
        
        with h5py.File(h5_filepath, 'r') as f:
            x_data = f['X_data'][:] 
            
            if 'SNR' in f:
                self.snr_data = f['SNR'][:]
            else:
                self.snr_data = np.zeros((x_data.shape[0], 1))
                
            # Simple Normalization
            flat = x_data.reshape(x_data.shape[0], -1)
            mins = flat.min(axis=1, keepdims=True)
            maxs = flat.max(axis=1, keepdims=True)
            rngs = maxs - mins
            rngs[rngs == 0] = 1.0
            mins = mins.reshape(-1, 1, 1)
            rngs = rngs.reshape(-1, 1, 1)
            x_data = (x_data - mins) / rngs
            self.data = torch.from_numpy(x_data).float().unsqueeze(1)

    def __len__(self):
        return len(self.sequence_starts)

    def __getitem__(self, idx):
        start = int(self.sequence_starts[idx])
        end = start + self.lookback_window
        seq = self.data[start:end] 
        label = self.labels[idx]
        
        snr_val = self.snr_data[start][0]
        
        return seq, torch.tensor(label, dtype=torch.long), torch.tensor(snr_val, dtype=torch.float)

def main():
    print("--- Testing Gold Code Transformer on UNSEEN SEEDS ---")
    
    raw_h5 = os.path.join('data', 'synthetic', 'classification_dataset_stft_random_deg10_python.h5')
    prep_h5 = os.path.join('data', 'synthetic', 'prepared_prediction_sequences_python.h5')
    
    if not os.path.exists(MODEL_PATH):
        print(f"Model file not found at {MODEL_PATH}")
        return

    if not os.path.exists(prep_h5):
        print(f"Prediction sequences not found at {prep_h5}")
        return

    # Load Test Meta (For Gold Code script, 'val' was basically 'unseen seeds'. 
    # But usually 'test' is even better. Let's check available keys.)
    with h5py.File(prep_h5, 'r') as f:
        # Looking for 'sequence_starts_test'
        if 'sequence_starts_test' in f:
            starts_test = f['sequence_starts_test'][:]
            y_test = f['YTest'][:]
            print("Using Test Set.")
        else:
            # Fallback to Validation if Test Not prepared (Old script logic)
            print("Test Set not found, using Validation Set instead.")
            starts_test = f['sequence_starts_validation'][:]
            y_test = f['YValidation'][:]
        
    print(f"Test Set Size: {len(starts_test)}")
    
    test_ds = RAMDataset(raw_h5, starts_test, y_test, LOOKBACK_WINDOW)
    test_loader = DataLoader(test_ds, batch_size=BATCH_SIZE, shuffle=False)
    
    # Load Model
    model = CNNTransformer(8).to(DEVICE)
    model.load_state_dict(torch.load(MODEL_PATH, map_location=DEVICE))
    model.eval()
    
    correct = 0
    total = 0
    all_preds = []
    all_labels = []
    
    # SNR Tracking
    snr_correct = {}
    snr_total = {}
    
    print("Running Inference...")
    with torch.no_grad():
        for X, y, snrs in test_loader:
            X, y, snrs = X.to(DEVICE), y.to(DEVICE), snrs.to(DEVICE)
            y = y.view(-1)
            
            out = model(X)
            pred = out.argmax(dim=1)
            
            correct += (pred == y).sum().item()
            total += y.size(0)
            
            # Bucketing
            pred_np = pred.cpu().numpy()
            y_np = y.cpu().numpy()
            snr_np = snrs.cpu().numpy()
            
            for i in range(len(y_np)):
                s = float(snr_np[i])
                if s not in snr_total:
                    snr_total[s] = 0
                    snr_correct[s] = 0
                snr_total[s] += 1
                if pred_np[i] == y_np[i]:
                    snr_correct[s] += 1
            
            all_preds.extend(pred_np)
            all_labels.extend(y_np)
            
    acc = 100 * correct / total
    print(f"Final TEST Accuracy: {acc:.2f}%")
    
    # Results per SNR
    sorted_snrs = sorted(snr_total.keys())
    snr_accs = []
    print("\n--- Accuracy per SNR ---")
    for s in sorted_snrs:
        ac = 100 * snr_correct[s] / snr_total[s]
        snr_accs.append(ac)
        print(f"SNR {s:5.1f} dB: {ac:5.1f}% ({snr_correct[s]}/{snr_total[s]})")
        
    # Plot SNR Curve
    plt.figure(figsize=(8, 6))
    plt.plot(sorted_snrs, snr_accs, marker='o', linewidth=2, color='orange')
    plt.title(f'Gold Code Prediction Accuracy vs SNR (Overall: {acc:.1f}%)')
    plt.xlabel('SNR (dB)')
    plt.ylabel('Accuracy (%)')
    plt.grid(True)
    plt.ylim(0, 100)
    plt.savefig('results/gold_accuracy_vs_snr.png')
    print("Saved results/gold_accuracy_vs_snr.png")
    
    # Confusion Matrix
    cm = confusion_matrix(all_labels, all_preds)
    plt.figure(figsize=(8, 6))
    sns.heatmap(cm, annot=True, fmt='d', cmap='Oranges')
    plt.title(f'Confusion Matrix (Acc: {acc:.1f}%)')
    plt.xlabel('Predicted')
    plt.ylabel('True')
    plt.savefig('results/test_gold_confusion.png')
    print("Saved results/test_gold_confusion.png")

if __name__ == "__main__":
    main()
