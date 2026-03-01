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
EPOCHS = 150
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
            
            # Load stats
            global_min = f['global_min'][0] if 'global_min' in f else -120.0
            global_max = f['global_max'][0] if 'global_max' in f else 0.0
            
            # Normalize Globally (Must match train_cnn.py EXACTLY!)
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
        label = self.labels[idx]
        
        return seq, torch.tensor(label, dtype=torch.long)

# --- 2. Positional Encoding & Periodic Activations ---
class PositionalEncoding(nn.Module):
    def __init__(self, d_model, max_len=5000):
        super(PositionalEncoding, self).__init__()
        # We don't want this in state dict so use register_buffer
        pe = torch.zeros(max_len, d_model)
        position = torch.arange(0, max_len, dtype=torch.float).unsqueeze(1)
        div_term = torch.exp(torch.arange(0, d_model, 2).float() * (-math.log(10000.0) / d_model))
        pe[:, 0::2] = torch.sin(position * div_term)
        pe[:, 1::2] = torch.cos(position * div_term)
        pe = pe.unsqueeze(0).transpose(0, 1) # -> (max_len, 1, d_model)
        self.register_buffer('pe', pe)

    def forward(self, x):
        # x shape should be (SeqLen, Batch, d_model)
        return x + self.pe[:x.size(0), :]


# THE SECRET SAUCE: Cosine Activation for Parity
class CosineActivation(nn.Module):
    def forward(self, x):
        # Using cos(pi * x) perfectly translates even/odd outputs into +1 / -1
        # This breaks the linear constraint of the ReLU!
        return torch.cos(math.pi * x)

# Custom Transformer Encoder Layer utilizing the Cosine Activation
class PeriodicTransformerEncoderLayer(nn.TransformerEncoderLayer):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        # Override the standard ReLU/GELU activation with our Periodic Cosine!
        self.activation = CosineActivation()


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
        
        # Project clean frequency probabilities (8 classes) to transformer hidden size
        self.prob_proj = nn.Linear(num_classes, d_model)
        
        # --- AUTOREGRESSIVE CO-T SUPERVISION ---
        self.fc = nn.Linear(d_model, 1)
        
        # Project Teacher Forcing bits
        self.bit_proj = nn.Linear(1, d_model)
        
        self.pos_encoder = PositionalEncoding(d_model)
        
        # --- THE FIX ---
        # Initialize our custom PERIODIC Transformer Layer instead of standard!
        encoder_layers = PeriodicTransformerEncoderLayer(d_model=d_model, nhead=nhead, dim_feedforward=128, dropout=0.0)
        self.transformer_encoder = nn.TransformerEncoder(encoder_layers, num_layers=num_layers)

    def forward(self, x, intermediate_bits=None):
        """
        x: (Batch, Seq, Channels, Height, Width) -> Standard STFT Hops
        intermediate_bits: (Batch, Num_Bits, 1) -> The Teacher Forcing bits.
        """
        b, s, c, h, w = x.size()
        c_in = x.view(b * s, c, h, w)
        
        c_out = self.cnn(c_in)
        c_flat = c_out.view(c_out.size(0), -1)
        
        logits = self.cnn_fc(c_flat)
        probs = torch.softmax(logits, dim=-1)
        
        src_spatial = self.prob_proj(probs).view(b, s, -1)
        
        if intermediate_bits is not None:
             src_bits = self.bit_proj(intermediate_bits)
             src = torch.cat([src_spatial, src_bits], dim=1)
             src = src.permute(1, 0, 2)
             src = self.pos_encoder(src)
             
             seq_len = src.size(0)
             mask = nn.Transformer.generate_square_subsequent_mask(seq_len).to(src.device)
             
             out = self.transformer_encoder(src, mask=mask)
             
             thought_outputs = out[s-1 : s-1 + 3, :, :]
             logits_out = self.fc(thought_outputs).squeeze(-1).permute(1, 0)
             return logits_out
             
        else:
            generated_bits = []
            current_bit_token = torch.zeros(b, 1, 1).to(x.device)
            
            for step in range(3):
                src_bit = self.bit_proj(current_bit_token)
                
                if step == 0:
                    src = src_spatial
                else:
                    src = torch.cat([src_spatial, src_bit], dim=1)
                
                src = src.permute(1, 0, 2)
                src = self.pos_encoder(src)
                
                seq_len = src.size(0)
                mask = nn.Transformer.generate_square_subsequent_mask(seq_len).to(src.device)
                
                out = self.transformer_encoder(src, mask=mask)
                
                next_bit_logit = self.fc(out[s-1 + step, :, :])
                next_bit = (next_bit_logit > 0).float()
                
                generated_bits.append(next_bit_logit)
                
                if step == 0:
                    current_bit_token = next_bit.unsqueeze(1)
                elif step < 2:
                     current_bit_token = torch.cat([current_bit_token, next_bit.unsqueeze(1)], dim=1)
            
            return torch.cat(generated_bits, dim=-1)

def main():
    print("--- Training PERIODIC CNN-Transformer on GOLD CODES ---")
    
    raw_h5 = os.path.join('data', 'synthetic', 'classification_dataset_stft_random_deg10_python.h5')
    prep_h5 = os.path.join('data', 'synthetic', 'prepared_prediction_sequences_python.h5')
    
    if not os.path.exists(raw_h5) or not os.path.exists(prep_h5):
        print(f"Data files missing!\nRaw: {raw_h5}\nPrep: {prep_h5}")
        return
        
    with h5py.File(prep_h5, 'r') as f:
        starts_train = f['sequence_starts_train'][:]
        y_train = f['YTrain'][:]
        starts_val = f['sequence_starts_validation'][:]
        y_val = f['YValidation'][:]
        
    print(f"Meta Loaded. Train: {len(starts_train)}, Val: {len(starts_val)}")
    
    train_ds = RAMDataset(raw_h5, starts_train, y_train, LOOKBACK_WINDOW)
    val_ds = RAMDataset(raw_h5, starts_val, y_val, LOOKBACK_WINDOW)
    
    train_loader = DataLoader(train_ds, batch_size=BATCH_SIZE, shuffle=True)
    val_loader = DataLoader(val_ds, batch_size=BATCH_SIZE, shuffle=False)
    
    model = CNNTransformer(8).to(DEVICE)
    
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
            elif k[0].isdigit():
                new_cnn_dict[k] = v
            elif k in ['weight', 'bias']:
                new_fc_dict[k] = v
                
        try:
            model.cnn.load_state_dict(new_cnn_dict, strict=True)
            if 'weight' in new_fc_dict and new_fc_dict['weight'].shape == model.cnn_fc.weight.shape:
                model.cnn_fc.weight.data = new_fc_dict['weight']
                model.cnn_fc.bias.data = new_fc_dict['bias']
            print("Loaded CNN Backbone + Classification Head successfully.")
        except RuntimeError as e:
            print(f"WARNING: Pretrained weights failed to load (Mismatch)! Error: {e}")
            print("Continuing training from scratch.")
    else:
        print("WARNING: Pretrained weights not found! Training from scratch.")
        
    for param in model.cnn.parameters():
        param.requires_grad = False
    for param in model.cnn_fc.parameters():
        param.requires_grad = False
    
    trainable_params = filter(lambda p: p.requires_grad, model.parameters())
    
    optimizer = optim.AdamW(trainable_params, lr=LEARNING_RATE)
    scheduler = optim.lr_scheduler.CosineAnnealingLR(optimizer, T_max=EPOCHS)
    criterion = nn.BCEWithLogitsLoss()
    
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
            
            y_bits = torch.zeros(y.size(0), 3).to(DEVICE)
            for i in range(3):
                y_bits[:, i] = (y >> (2 - i)) & 1
            
            teacher_bits = y_bits[:, :2].unsqueeze(-1)
            
            optimizer.zero_grad()
            out = model(X, intermediate_bits=teacher_bits)
            
            loss = criterion(out, y_bits)
            loss.backward()
            optimizer.step()
            
            total_loss += loss.item()
            
            pred_bits = (out > 0).float()
            correct += (pred_bits == y_bits).sum().item()
            total += y.size(0) * 3 
            
            correct_seq += (pred_bits == y_bits).all(dim=1).sum().item()
            total_seq += y.size(0)
            
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
                
                out = model(X, intermediate_bits=None)
                v_loss_batch = criterion(out, y_bits)
                val_loss += v_loss_batch.item()
                
                pred_bits = (out > 0).float()
                v_corr += (pred_bits == y_bits).sum().item()
                v_tot += y.size(0) * 3
                
                v_seq_corr += (pred_bits == y_bits).all(dim=1).sum().item()
                v_seq_tot += y.size(0)
        
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
            output_model = os.path.join('models', 'gold_code_periodic_transformer_best.pth')
            # create dir if not exists
            os.makedirs(os.path.dirname(output_model), exist_ok=True)
            torch.save(model.state_dict(), output_model)
            print("  --> Best Model Saved")
            
    plt.figure(figsize=(12, 5))
    
    plt.subplot(1, 2, 1)
    plt.plot(history['train_loss'], label='Train Loss')
    plt.plot(history['val_loss'], label='Val Loss')
    plt.title('Periodic Gold Code Loss History')
    plt.xlabel('Epoch')
    plt.legend()
    
    plt.subplot(1, 2, 2)
    plt.plot(history['train_acc'], label='Train Acc')
    plt.plot(history['val_acc'], label='Val Acc')
    plt.title('Periodic Gold Code Accuracy History')
    plt.xlabel('Epoch')
    plt.legend()
    
    plt.tight_layout()
    # create dir if not exists
    os.makedirs('results', exist_ok=True)
    plt.savefig('results/training_history_periodic_gold.png')
    print("Training Complete. History plot saved to results/training_history_periodic_gold.png")

if __name__ == "__main__":
    main()
