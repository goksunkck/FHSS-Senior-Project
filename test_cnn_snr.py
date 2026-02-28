
import torch
import torch.nn as nn
from torch.utils.data import Dataset, DataLoader
import h5py
import numpy as np
import os
import matplotlib.pyplot as plt

# --- Same Config ---
DEVICE = torch.device("cuda" if torch.cuda.is_available() else "cpu")
BATCH_SIZE = 256

class SNRDataset(Dataset):
    def __init__(self, h5_path):
        print(f"Loading {h5_path}...")
        with h5py.File(h5_path, 'r') as f:
            self.X = f['X_data'][:]
            self.Y = f['Y_data'][:]
            if 'SNR' in f:
                self.SNR = f['SNR'][:]
            else:
                print("Warning: No SNR data found, assuming 0dB")
                self.SNR = np.zeros_like(self.Y)
            
            # Simple Norm
            flat = self.X.reshape(self.X.shape[0], -1)
            mins = flat.min(axis=1, keepdims=True)
            maxs = flat.max(axis=1, keepdims=True)
            rngs = maxs - mins
            rngs[rngs == 0] = 1.0
            mins = mins.reshape(-1, 1, 1)
            rngs = rngs.reshape(-1, 1, 1)
            self.X = (self.X - mins) / rngs
            
            self.X = torch.from_numpy(self.X).float()
            self.Y = torch.from_numpy(self.Y).long().squeeze()
            self.SNR = torch.from_numpy(self.SNR).float().squeeze()

    def __len__(self):
        return self.X.shape[0]

    def __getitem__(self, idx):
        return self.X[idx].unsqueeze(0), self.Y[idx], self.SNR[idx]

# --- Re-Define Classes to Load Model ---
class PositionalEncoding(nn.Module):
    def __init__(self, d_model, max_len=5000):
        super(PositionalEncoding, self).__init__()
        pe = torch.zeros(max_len, d_model)
        position = torch.arange(0, max_len, dtype=torch.float).unsqueeze(1)
        div_term = torch.exp(torch.arange(0, d_model, 2).float() * (-np.log(10000.0) / d_model))
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
        
        self.embedding = nn.Embedding(num_classes, d_model)
        self.pos_encoder = PositionalEncoding(d_model)
        
        encoder_layers = nn.TransformerEncoderLayer(d_model=d_model, nhead=nhead, dim_feedforward=128, dropout=0.0)
        self.transformer_encoder = nn.TransformerEncoder(encoder_layers, num_layers=num_layers)
        
        self.fc = nn.Linear(d_model, num_classes)

    def forward(self, x):
        # We only need the vision part for this test
        b, s, c, h, w = x.size() # But here inputs are (N, 1, 1024, 3) from SingleFrame loader
        # actually, the main forward expects sequence.
        # Let's define a helper to just run vision
        pass

    def vision_forward(self, x):
        # x: (N, 1, 1024, 3)
        c_out = self.cnn(x)
        c_flat = c_out.view(c_out.size(0), -1)
        logits = self.cnn_fc(c_flat)
        return logits

def main():
    print("--- Testing CNN Consistency vs SNR ---")
    h5_path = os.path.join('data', 'synthetic', 'm_sequence_dataset.h5')
    # We load the TRANSFORMER model because it has the exact CNN + Random Head state used during testing
    model_path = os.path.join('models', 'm_sequence_transformer_best.pth')
    
    if not os.path.exists(h5_path):
        print("Dataset missing.")
        return
        
    ds = SNRDataset(h5_path)
    # Use a subset 
    indices = torch.randperm(len(ds))[:50000] 
    subset = torch.utils.data.Subset(ds, indices)
    loader = DataLoader(subset, batch_size=BATCH_SIZE, shuffle=False)
    
    model = CNNTransformer(8).to(DEVICE)
    if os.path.exists(model_path):
        print(f"Loading weights from {model_path}")
        model.load_state_dict(torch.load(model_path, map_location=DEVICE))
    else:
        print("Model missing!")
        return

    model.eval()
    
    # We need to map: SNR -> GroundTruth -> [Predictions List]
    # e.g. results[-10.0][5] = [3, 3, 4, 1, 3]
    results = {} 
    
    print("Evaluating...")
    with torch.no_grad():
        for X, y, snr in loader:
            X, y, snr = X.to(DEVICE), y.to(DEVICE), snr.to(DEVICE)
            
            # Use vision forward
            logits = model.vision_forward(X)
            pred = logits.argmax(dim=1)
            
            y_np = y.cpu().numpy()
            snr_np = snr.cpu().numpy()
            pred_np = pred.cpu().numpy()
            
            for i in range(len(y_np)):
                s = float(snr_np[i])
                lbl = int(y_np[i])
                p = int(pred_np[i])
                
                if s not in results: results[s] = {}
                if lbl not in results[s]: results[s][lbl] = []
                results[s][lbl].append(p)

    # Analyze Consistency
    snrs = sorted(results.keys())
    consistency_scores = []
    
    print("\n--- CNN Consistency (Effective Accuracy) vs SNR ---")
    for s in snrs:
        total_correct_mappings = 0
        total_samples = 0
        
        # For each ground truth label, find the 'Mode' (Most frequent prediction)
        # If the model is consistent, it should always predict the Mode.
        for lbl in range(8):
            if lbl not in results[s]: continue
            preds = results[s][lbl]
            if len(preds) == 0: continue
            
            # Find Mode
            counts = np.bincount(preds, minlength=8)
            mode = np.argmax(counts)
            
            # Count how many times it predicted the Mode
            total_correct_mappings += counts[mode]
            total_samples += len(preds)
            
        consistency = 100 * total_correct_mappings / total_samples
        consistency_scores.append(consistency)
        print(f"SNR {s:5.1f} dB: {consistency:5.1f}%")
        
    plt.figure(figsize=(8,6))
    plt.plot(snrs, consistency_scores, marker='o', linewidth=2, color='green')
    plt.title("CNN Effective 'Consistency' vs SNR\n(Correcting for Random Initialized Head)")
    plt.xlabel("SNR (dB)")
    plt.ylabel("Consistency (%)")
    plt.grid(True)
    plt.ylim(0, 105)
    plt.savefig('results/cnn_consistency_vs_snr.png')
    print("Saved results/cnn_consistency_vs_snr.png")

if __name__ == "__main__":
    main()
