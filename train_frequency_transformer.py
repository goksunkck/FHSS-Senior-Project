
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
        
        # Project clean frequency probabilities (8 classes) to transformer hidden size
        self.prob_proj = nn.Linear(num_classes, d_model)
        
        # --- AUTOREGRESSIVE CO-T SUPERVISION ---
        # The paper states: "The model must generate CoT chains end-to-end, parity can be learned 
        # efficiently if augmented data is employed to internally verify the soundness of 
        # intermediate steps" AND "intermediate parities are incorporated into the loss".
        # To accomplish "Teacher Forcing / Intermediate Supervision", we add a recurrent
        # output sequence. We output a dimension of 1 (a single bit) iteratively.
        self.fc = nn.Linear(d_model, 1)
        
        # We also need to project the previously generated bits back into the transformer 
        # to act as the "Teacher Forcing" tokens!
        self.bit_proj = nn.Linear(1, d_model)
        
        # We must NOT use batch_first=True natively because we need causal masking 
        # so the model can't cheat and look ahead. Standard (Seq, Batch, Embed) is easier.
        self.pos_encoder = PositionalEncoding(d_model)
        encoder_layers = nn.TransformerEncoderLayer(d_model=d_model, nhead=nhead, dim_feedforward=128, dropout=0.0)
        self.transformer_encoder = nn.TransformerEncoder(encoder_layers, num_layers=num_layers)

    def forward(self, x, intermediate_bits=None):
        """
        x: (Batch, Seq, Channels, Height, Width) -> Standard STFT Hops
        intermediate_bits: (Batch, Num_Bits, 1) -> The Teacher Forcing bits.
                           If None, we do inference autoregressively!
        """
        b, s, c, h, w = x.size()
        c_in = x.view(b * s, c, h, w)
        
        # Get continuous features from the CNN
        c_out = self.cnn(c_in)
        c_flat = c_out.view(c_out.size(0), -1)
        
        # Pass through the pre-trained classification head
        logits = self.cnn_fc(c_flat)
        
        # Convert to a soft probability distribution (differentiable and clean!)
        probs = torch.softmax(logits, dim=-1)
        
        # Project clean probabilities to d_model space: (Batch * SeqLen, 8) -> (Batch * SeqLen, d_model)
        # Reshape to (Batch, SeqLen, d_model)
        src_spatial = self.prob_proj(probs).view(b, s, -1)
        
        # --- TEACHER FORCING LOOP (CHAIN OF THOUGHT) ---
        # Instead of doing 1 forward pass, we create a sequence of thoughts!
        # [Spatial Hop 1] ... [Spatial Hop 12] + [Bit 1 Thought] + [Bit 2 Thought] -> [Predicts Bit 3]
        
        if intermediate_bits is not None:
             # Training Mode (Teacher Forcing)
             # We inject the true intermediate bits as context!
             
             # (Batch, Num_Bits, d_model)
             src_bits = self.bit_proj(intermediate_bits)
             
             # Concatenate Temporal Sequence + Thought Sequence
             # (Batch, SeqLen + Num_Bits, d_model)
             src = torch.cat([src_spatial, src_bits], dim=1)
             
             # Convert to Transformer format: (SeqLen, Batch, d_model)
             src = src.permute(1, 0, 2)
             src = self.pos_encoder(src)
             
             # Causal Mask to prevent cheating (predicting bit 2 while looking at bit 3)
             seq_len = src.size(0)
             mask = nn.Transformer.generate_square_subsequent_mask(seq_len).to(src.device)
             
             out = self.transformer_encoder(src, mask=mask)
             
             # We only care about the outputs corresponding to the "Thought" tokens
             # We want to predict Bit 1, Bit 2, Bit 3.
             # Output at index 'Spatial_Sequence_Length - 1' predicts Bit 1
             # Output at index 'Spatial_Sequence_Length' predicts Bit 2...
             
             thought_outputs = out[s-1 : s-1 + intermediate_bits.size(1), :, :]
             
             # Output shape: (Num_Bits, Batch, 1) -> (Batch, Num_Bits)
             logits_out = self.fc(thought_outputs).squeeze(-1).permute(1, 0)
             return logits_out
             
        else:
            # Inference Mode (End-to-End CoT Generation)
            generated_bits = []
            
            # Start with an empty thought vector (e.g. 0) or just the spatial sequence
            # For simplicity, we just use 0 as the "Start" bit token internally.
            current_bit_token = torch.zeros(b, 1, 1).to(x.device)
            
            for step in range(3):
                # (Batch, 1, d_model)
                src_bit = self.bit_proj(current_bit_token)
                
                # Combine what we have so far
                if step == 0:
                    src = src_spatial
                else:
                    # Append the thought sequence
                    # e.g [Spatial...] + [Bit 1]
                    src = torch.cat([src_spatial, src_bit], dim=1)
                
                src = src.permute(1, 0, 2)
                src = self.pos_encoder(src)
                
                # We don't technically need a mask during step-by-step inference, 
                # but standard practice is to use it.
                seq_len = src.size(0)
                mask = nn.Transformer.generate_square_subsequent_mask(seq_len).to(src.device)
                
                out = self.transformer_encoder(src, mask=mask)
                
                # The very last token in the sequence predicts the next piece of parity math
                last_token = out[-1, :, :]
                
                # Predict 1 Parity Bit
                next_bit_logit = self.fc(last_token) 
                
                # Hard decision for the next teacher forcing loop!
                next_bit = (next_bit_logit > 0).float()
                
                generated_bits.append(next_bit_logit)
                
                # Append this hard decision to our thought stream for the next iteration
                if step == 0:
                    current_bit_token = next_bit.unsqueeze(1)
                else:
                     current_bit_token = torch.cat([current_bit_token, next_bit.unsqueeze(1)], dim=1)
            
            # (Batch, 3)
            return torch.cat(generated_bits, dim=-1)

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
        
        print(f"DEBUG: Keys in loaded state_dict: {list(d.keys())}")
        
        for k, v in d.items():
            if k.startswith('cnn.'):
                name = k[4:] 
                new_cnn_dict[name] = v
            elif k.startswith('fc.'):
                name = k[3:]
                new_fc_dict[name] = v
            elif k[0].isdigit():
                new_cnn_dict[k] = v
            elif k == 'weight' or k == 'bias':
                new_fc_dict[k] = v
                
        try:
            model.cnn.load_state_dict(new_cnn_dict, strict=True)
            # Make sure we only load fc weights if shape matches
            if 'weight' in new_fc_dict and new_fc_dict['weight'].shape == model.cnn_fc.weight.shape:
                model.cnn_fc.weight.data = new_fc_dict['weight']
                model.cnn_fc.bias.data = new_fc_dict['bias']
            print("Loaded CNN Backbone + Classification Head successfully.")
        except RuntimeError as e:
            print(f"WARNING: Pretrained weights failed to load (Mismatch)! Error: {e}")
            print("Continuing training from scratch.")
    else:
        print("WARNING: Pretrained weights not found! Training from scratch.")
        
    # FREEZE the CNN Backbone AND the classification head!
    for param in model.cnn.parameters():
        param.requires_grad = False
    for param in model.cnn_fc.parameters():
        param.requires_grad = False
    
    # Only optimize the remaining layers (Transformer + Projections + Heads)
    trainable_params = filter(lambda p: p.requires_grad, model.parameters())
    
    optimizer = optim.AdamW(trainable_params, lr=LEARNING_RATE)
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
            
            # Create the Teacher Forcing Input Sequence
            # We want to use a "Start" token (Bit=0), then True Bit 1, then True Bit 2.
            # We DONT input True Bit 3, because the model is predicting Bit 3!
            
            # Shape (Batch, 3)
            # Intermediate teacher bits: [0, y1, y2]
            teacher_bits = torch.zeros(y.size(0), 3).to(DEVICE)
            teacher_bits[:, 1] = y_bits[:, 0]
            teacher_bits[:, 2] = y_bits[:, 1]
            teacher_bits = teacher_bits.unsqueeze(-1) # (Batch, 3, 1)
            
            optimizer.zero_grad()
            out = model(X, intermediate_bits=teacher_bits)
            
            # --- CHAIN OF THOUGHT LOSS ---
            # Paper: Parity can be learned efficiently when intermediate parities 
            # are incorporated into the loss function ("Teacher Forcing").
            # Here, we penalize the model for getting ANY of the 3 intermediate LFSR bits wrong
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
                
                # INFERENCE MODE: No teacher forcing! The model must generate CoT chains end-to-end
                out = model(X, intermediate_bits=None)
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
