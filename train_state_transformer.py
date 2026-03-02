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
EPOCHS = 350
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
        
        # --- PRECOMPUTE LFSR STATES FOR O(1) LOOKUP ---
        # Instead of computing the LFSR step inside __getitem__ which is incredibly slow,
        # we precompute the 256 states for all 1024 unique Gold Code seeds!
        print("Precomputing LFSR internal states for all seeds to speed up training...")
        self.state_cache = {}
        for seed_int in range(1, 1024):  # 1023 max seeds
            pn_initial1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 1]
            seed_bin_str = format(seed_int, '010b')
            pn_initial2 = [int(b) for b in seed_bin_str]
            
            lfsr1 = SimpleLFSR([10, 3, 0], pn_initial1)
            lfsr2 = SimpleLFSR([10, 7, 0], pn_initial2)
            
            states_for_seed = []
            for hop in range(256):
                # Save state BEFORE advancing for this hop
                state1 = lfsr1.get_state_bits()
                state2 = lfsr2.get_state_bits()
                
                # Convert to -1, +1 immediately for our math
                bipolar_state = [(b * 2.0 - 1.0) for b in (state1 + state2)]
                states_for_seed.append(torch.tensor(bipolar_state, dtype=torch.float))
                
                # Advance 3 bits for the next hop
                lfsr1.step()
                lfsr2.step()
                lfsr1.step()
                lfsr2.step()
                lfsr1.step()
                lfsr2.step()
                
            self.state_cache[seed_int] = states_for_seed
            
        # For seed 0, which might technically exist but is degenerate
        self.state_cache[0] = [torch.zeros(20, dtype=torch.float) for _ in range(256)]
        
        print("LFSR Cache built successfully!")

    def __len__(self):
        return len(self.sequence_starts)

    def __getitem__(self, idx):
        start = int(self.sequence_starts[idx])
        end = start + self.lookback_window
        seq = self.data[start:end] 
        
        # --- O(1) LOOKUP TRUE INTERNAL STATE ---
        target_global_idx = end # The hop we are trying to predict
        local_target_hop = target_global_idx % 256 # Assuming num_hops = 256
        seed_int = int(self.seed_vals[target_global_idx][0])
        
        # Fetch the exact 20-bit internal memory state at this exact moment!
        label = self.state_cache[seed_int][local_target_hop]
        
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
        self.bit_proj = nn.Linear(1, d_model)
        
        self.pos_encoder = PositionalEncoding(d_model)
        
        # Periodic Transformer
        encoder_layers = PeriodicTransformerEncoderLayer(d_model=d_model, nhead=nhead, dim_feedforward=128, dropout=0.0)
        self.transformer_encoder = nn.TransformerEncoder(encoder_layers, num_layers=num_layers)
        
        # Step-by-step logic predicts 1 bit at a time
        self.fc = nn.Linear(d_model, 1)

    def generate_square_subsequent_mask(self, sz):
        mask = (torch.triu(torch.ones(sz, sz)) == 1).transpose(0, 1)
        mask = mask.float().masked_fill(mask == 0, float('-inf')).masked_fill(mask == 1, float(0.0))
        return mask

    def forward(self, x, y_state=None):
        import torch.nn.functional as F
        b, s, c, h, w = x.size()
        c_in = x.view(b * s, c, h, w)
        
        # Freeze CNN parts and force EXACT discrete boolean values
        with torch.no_grad():
            c_out = self.cnn(c_in)
            c_flat = c_out.view(c_out.size(0), -1)
            logits = self.cnn_fc(c_flat)
            # Emulate pure logical inputs via argmax instead of fuzzy continuous softmax probabilities
            preds = torch.argmax(logits, dim=-1)
            preds_onehot = F.one_hot(preds, num_classes=8).float()
            
        vis_emb = self.prob_proj(preds_onehot).view(b, s, -1)
        
        if y_state is not None:
            # --- TEACHER FORCING MULTI-STEP (Chain of Thought training) ---
            teacher_bits = y_state[:, :-1].unsqueeze(-1)
            state_emb = self.bit_proj(teacher_bits)
            
            src = torch.cat([vis_emb, state_emb], dim=1)
            seq_len = src.size(1)
            
            src = src.permute(1, 0, 2)
            src = self.pos_encoder(src)
            
            mask = self.generate_square_subsequent_mask(seq_len).to(src.device)
            out = self.transformer_encoder(src, mask=mask)
            
            # Extract exactly the 20 output elements corresponding to the predictions
            state_preds = self.fc(out[s-1:, :, :])
            return state_preds.squeeze(-1).permute(1, 0)
        else:
            # --- AUTOREGRESSIVE INFERENCE ---
            curr_src = vis_emb
            generated_bits = []
            
            for step in range(20):
                src_perm = curr_src.permute(1, 0, 2)
                src_encoded = self.pos_encoder(src_perm)
                
                seq_len = src_encoded.size(0)
                mask = self.generate_square_subsequent_mask(seq_len).to(curr_src.device)
                
                out = self.transformer_encoder(src_encoded, mask=mask)
                
                next_bit_logit = self.fc(out[-1, :, :])
                generated_bits.append(next_bit_logit)
                
                if step < 19:
                    next_bit = (next_bit_logit > 0).float() * 2.0 - 1.0
                    next_bit_emb = self.bit_proj(next_bit).unsqueeze(1)
                    curr_src = torch.cat([curr_src, next_bit_emb], dim=1)
                    
            return torch.cat(generated_bits, dim=-1)

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
    
    # Set layers and heads to 1 for exact mathematical match to proofs
    model = StateTransformer(num_classes=8, d_model=64, nhead=4, num_layers=2).to(DEVICE)
    
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
    
    # MSE loss is perfect for targets mapped to -1 and 1
    criterion = nn.MSELoss()
    
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
            out = model(X, y_state=y_state) # Teacher Forcing properly applied!
            
            loss = criterion(out, y_state)
            loss.backward()
            optimizer.step()
            
            total_loss += loss.item()
            
            # Map predictions back to (-1, 1) space
            pred_bits = (out > 0).float() * 2.0 - 1.0
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
                
                # Do not pass y_state during validation to force True Autoregressive Inference
                out = model(X)
                v_loss_batch = criterion(out, y_state)
                val_loss += v_loss_batch.item()
                
                # Map evaluation predictions to {-1, 1} space
                pred_bits = (out > 0).float() * 2.0 - 1.0
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
