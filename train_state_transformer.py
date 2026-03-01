import torch
import torch.nn as nn
import torch.optim as optim
from torch.utils.data import Dataset, DataLoader
import numpy as np
import h5py
import os
import math
import matplotlib.pyplot as plt
import sys

# --- Configuration ---
BATCH_SIZE = 128
EPOCHS = 150
LEARNING_RATE = 1e-3
LOOKBACK_WINDOW = 12
DEVICE = torch.device("cuda" if torch.cuda.is_available() else "cpu")

# --- LFSR Helper for True State Generation ---
class SimpleLFSR:
    def __init__(self, poly, initial_state):
        self.degree = poly[0]
        self.taps = [p for p in poly if 0 < p < self.degree]
        self.state = 0
        for b in initial_state:
            self.state = (self.state << 1) | int(b)
            
    def step(self):
        out_bit = (self.state >> (self.degree - 1)) & 1
        feedback = out_bit
        for t in self.taps:
            feedback ^= (self.state >> (t - 1)) & 1
        self.state = ((self.state << 1) & ((1 << self.degree) - 1)) | feedback
        return out_bit
        
    def get_state_bits(self):
        return [(self.state >> i) & 1 for i in range(self.degree)]

# --- 1. Dataset (RAM Loading with True State Labeling) ---
class RAMStateDataset(Dataset):
    def __init__(self, h5_filepath, sequence_starts, lookback_window):
        print(f"preload: Loading dataset from {h5_filepath} into RAM...")
        self.sequence_starts = sequence_starts
        self.lookback_window = lookback_window
        
        with h5py.File(h5_filepath, 'r') as f:
            x_data = f['X_data'][:] 
            self.seed_vals = f['Seed_Val'][:] # We need the seeds to recreate the internal state!
            
            # Load stats
            global_min = f['global_min'][0] if 'global_min' in f else -120.0
            global_max = f['global_max'][0] if 'global_max' in f else 0.0
            
            # Normalize Globally
            print(f"Normalizing globally... Min: {global_min}, Max: {global_max}")
            x_data = (x_data - global_min) / (global_max - global_min + 1e-6)
            self.data = torch.from_numpy(x_data).float().unsqueeze(1) # (N, 1, 1024, 3)
            
            print(f"Loaded {self.data.shape}")

    def __len__(self):
        return len(self.sequence_starts)

    def __getitem__(self, idx):
        start = int(self.sequence_starts[idx])
        end = start + self.lookback_window
        seq = self.data[start:end] 
        
        # --- GENERATING THE TRUE INTERNAL STATE ---
        target_global_idx = end # The hop we are trying to predict
        local_target_hop = target_global_idx % 256 # Assuming num_hops = 256
        
        seed_int = int(self.seed_vals[target_global_idx][0])
        
        # Init Gold Code LFSRs
        pn_initial1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 1]
        seed_bin_str = format(seed_int, '010b')
        pn_initial2 = [int(b) for b in seed_bin_str]
        
        lfsr1 = SimpleLFSR([10, 3, 0], pn_initial1)
        lfsr2 = SimpleLFSR([10, 7, 0], pn_initial2)
        
        # Advance the LFSR exactly up to the moment BEFORE this hop is generated
        # 3 bits generated per hop
        bits_to_advance = local_target_hop * 3
        for _ in range(bits_to_advance):
            lfsr1.step()
            lfsr2.step()
            
        # The true 20-bit internal memory state of the device at this exact moment!
        # This is the ultimate "Chain of Thought" teacher.
        state1 = lfsr1.get_state_bits()
        state2 = lfsr2.get_state_bits()
        full_20_bit_state = state1 + state2
        
        label = torch.tensor(full_20_bit_state, dtype=torch.float)
        
        return seq, label

# --- 2. Positional Encoding & Periodic Activations ---
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

class CosineActivation(nn.Module):
    def forward(self, x):
        return torch.cos(math.pi * x)

class PeriodicTransformerEncoderLayer(nn.TransformerEncoderLayer):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.activation = CosineActivation()

# --- 3. Model Architecture ---
class StateTransformer(nn.Module):
    def __init__(self, num_classes=8, d_model=64, nhead=4, num_layers=2):
        super(StateTransformer, self).__init__()
        
        # Vision
        self.cnn = nn.Sequential(
            nn.Conv2d(1, 32, kernel_size=(5, 1), stride=1, padding=(2, 0)),
            nn.InstanceNorm2d(32),
            nn.ReLU(),
            nn.MaxPool2d(kernel_size=(4, 1), stride=(4, 1)),
            nn.Conv2d(32, 64, kernel_size=(5, 1), stride=1, padding=(2, 0)),
            nn.InstanceNorm2d(64),
            nn.ReLU(),
            nn.MaxPool2d(kernel_size=(4, 1), stride=(4, 1)),
            nn.Conv2d(64, 32, kernel_size=(5, 1), stride=1, padding=(2, 0)),
            nn.InstanceNorm2d(32),
            nn.ReLU(),
            nn.Conv2d(32, 32, kernel_size=(5, 1), stride=1, padding=(2, 0)),
            nn.InstanceNorm2d(32),
            nn.ReLU()
        )
        
        self.flat_size = 32 * 64 * 3 
        self.cnn_fc = nn.Linear(self.flat_size, num_classes)
        self.prob_proj = nn.Linear(num_classes, d_model)
        
        self.pos_encoder = PositionalEncoding(d_model)
        
        # Periodic Transformer
        encoder_layers = PeriodicTransformerEncoderLayer(d_model=d_model, nhead=nhead, dim_feedforward=128, dropout=0.0)
        self.transformer_encoder = nn.TransformerEncoder(encoder_layers, num_layers=num_layers)
        
        # Output is NOT 3 bits anymore! It predicts the full 20 bit internal state of the device!
        self.fc = nn.Linear(d_model, 20)

    def forward(self, x):
        b, s, c, h, w = x.size()
        c_in = x.view(b * s, c, h, w)
        
        c_out = self.cnn(c_in)
        c_flat = c_out.view(c_out.size(0), -1)
        
        logits = self.cnn_fc(c_flat)
        probs = torch.softmax(logits, dim=-1)
        
        src = self.prob_proj(probs).view(b, s, -1)
        
        src = src.permute(1, 0, 2)
        src = self.pos_encoder(src)
        
        # Because we predict everything in 1 shot, we don't strictly need a causal mask for an autoregressive bit loop.
        # But we DO still mask the spatial tokens depending on standard transformer practices.
        # But actually, PyTorch transformer by default uses all tokens if no mask is provided.
        out = self.transformer_encoder(src)
        
        # Take the final spatial sequence token to predict the state
        last_token = out[-1, :, :]
        
        # Output shape: (Batch, 20)
        return self.fc(last_token)

def main():
    print("--- Training State-Tracking Periodic Transformer ---")
    
    raw_h5 = os.path.join('data', 'synthetic', 'classification_dataset_stft_random_deg10_python.h5')
    prep_h5 = os.path.join('data', 'synthetic', 'prepared_prediction_sequences_python.h5')
    
    if not os.path.exists(raw_h5) or not os.path.exists(prep_h5):
        print(f"Data files missing!\nRaw: {raw_h5}\nPrep: {prep_h5}")
        return
        
    with h5py.File(prep_h5, 'r') as f:
        starts_train = f['sequence_starts_train'][:]
        starts_val = f['sequence_starts_validation'][:]
        
    print(f"Meta Loaded. Train: {len(starts_train)}, Val: {len(starts_val)}")
    
    train_ds = RAMStateDataset(raw_h5, starts_train, LOOKBACK_WINDOW)
    val_ds = RAMStateDataset(raw_h5, starts_val, LOOKBACK_WINDOW)
    
    train_loader = DataLoader(train_ds, batch_size=BATCH_SIZE, shuffle=True)
    val_loader = DataLoader(val_ds, batch_size=BATCH_SIZE, shuffle=False)
    
    model = StateTransformer(8).to(DEVICE)
    
    w_path = os.path.join('models', 'cnn_weights.pth')
    if os.path.exists(w_path):
        print("Loading Pretrained CNN (Full)...")
        d = torch.load(w_path, map_location=DEVICE)
        new_cnn_dict = {}
        new_fc_dict = {}
        for k, v in d.items():
            if k.startswith('cnn.'): new_cnn_dict[k[4:]] = v
            elif k.startswith('fc.'): new_fc_dict[k[3:]] = v
            elif k[0].isdigit(): new_cnn_dict[k] = v
            elif k in ['weight', 'bias']: new_fc_dict[k] = v
                
        try:
            model.cnn.load_state_dict(new_cnn_dict, strict=True)
            if 'weight' in new_fc_dict and new_fc_dict['weight'].shape == model.cnn_fc.weight.shape:
                model.cnn_fc.weight.data = new_fc_dict['weight']
                model.cnn_fc.bias.data = new_fc_dict['bias']
            print("Loaded CNN Backbone + Classification Head successfully.")
        except RuntimeError as e:
            print(f"WARNING: Pretrained weights failed to load (Mismatch)! Error: {e}")
    else:
        print("WARNING: Pretrained weights not found! Training from scratch.")
        
    for param in model.cnn.parameters(): param.requires_grad = False
    for param in model.cnn_fc.parameters(): param.requires_grad = False
    
    trainable_params = filter(lambda p: p.requires_grad, model.parameters())
    optimizer = optim.AdamW(trainable_params, lr=LEARNING_RATE)
    scheduler = optim.lr_scheduler.CosineAnnealingLR(optimizer, T_max=EPOCHS)
    
    # BCE loss is perfect for 20 independent bits!
    criterion = nn.BCEWithLogitsLoss()
    
    history = {'train_loss': [], 'val_loss': [], 'train_seq_acc': [], 'val_seq_acc': [], 'val_bit_acc': []}
    best_val_loss = float('inf')
    
    for epoch in range(EPOCHS):
        model.train()
        total_loss = 0
        correct_state = 0
        total_seq = 0
        
        for X, y_state in train_loader:
            X, y_state = X.to(DEVICE), y_state.to(DEVICE)
            
            optimizer.zero_grad()
            out = model(X) # Predict 20 bits
            
            loss = criterion(out, y_state)
            loss.backward()
            optimizer.step()
            
            total_loss += loss.item()
            
            pred_bits = (out > 0).float()
            # How many sequences did we nail the ENTIRE 20-bit state exactly correct?
            correct_state += (pred_bits == y_state).all(dim=1).sum().item()
            total_seq += y_state.size(0)
            
        val_loss = 0
        model.eval()
        v_seq_corr = 0
        v_seq_tot = 0
        v_bit_corr = 0
        v_bit_tot = 0
        with torch.no_grad():
            for X, y_state in val_loader:
                X, y_state = X.to(DEVICE), y_state.to(DEVICE)
                
                out = model(X)
                v_loss_batch = criterion(out, y_state)
                val_loss += v_loss_batch.item()
                
                pred_bits = (out > 0).float()
                v_seq_corr += (pred_bits == y_state).all(dim=1).sum().item()
                v_seq_tot += y_state.size(0)
                
                v_bit_corr += (pred_bits == y_state).sum().item()
                v_bit_tot += y_state.size(0) * 20
        
        avg_train_loss = total_loss / len(train_loader)
        avg_val_loss = val_loss / len(val_loader)
        train_seq_acc = 100 * correct_state / total_seq
        val_seq_acc = 100 * v_seq_corr / v_seq_tot
        val_bit_acc = 100 * v_bit_corr / v_bit_tot
        
        history['train_loss'].append(avg_train_loss)
        history['val_loss'].append(avg_val_loss)
        history['train_seq_acc'].append(train_seq_acc)
        history['val_seq_acc'].append(val_seq_acc)
        history['val_bit_acc'].append(val_bit_acc)
        
        scheduler.step()
        
        print(f"Epoch {epoch+1}/{EPOCHS} | Train Loss: {avg_train_loss:.4f} | Val Loss: {avg_val_loss:.4f} | Val Bit Acc: {val_bit_acc:.1f}% | Exact 20-Bit State Acc: {val_seq_acc:.1f}%")
        
        if avg_val_loss < best_val_loss:
            best_val_loss = avg_val_loss
            output_model = os.path.join('models', 'state_tracking_transformer_best.pth')
            os.makedirs(os.path.dirname(output_model), exist_ok=True)
            torch.save(model.state_dict(), output_model)
            if epoch > 0: print("  --> Best Model Saved")

if __name__ == "__main__":
    main()
