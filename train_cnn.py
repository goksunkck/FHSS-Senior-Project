import torch
import torch.nn as nn
import torch.optim as optim
from torch.utils.data import Dataset, DataLoader
import h5py
import numpy as np
import os

# --- Configuration ---
BATCH_SIZE = 32
EPOCHS = 10
LEARNING_RATE = 1e-3
DEVICE = torch.device("cuda" if torch.cuda.is_available() else "cpu")

# --- Dataset ---
class SingleFrameDataset(Dataset):
    def __init__(self, h5_path):
        self.h5_path = h5_path
        print(f"Loading dataset into RAM: {h5_path}...")
        with h5py.File(self.h5_path, 'r') as f:
            # Load all data into memory
            self.X = f['X_data'][:] # (N, 1024, 3)
            self.Y = f['Y_data'][:] # (N, 1)
            
            # Load stats
            self.global_min = f['global_min'][0] if 'global_min' in f else -120.0
            self.global_max = f['global_max'][0] if 'global_max' in f else 0.0
            
        print(f"Loaded {self.X.shape[0]} samples. Normalizing...")
        # Normalize in-place to save memory if possible, or just on the fly
        # Validating normalization logic: 0-1 range
        self.X = (self.X - self.global_min) / (self.global_max - self.global_min + 1e-6)
        
        # Convert to tensor
        self.X = torch.from_numpy(self.X).float()
        self.Y = torch.from_numpy(self.Y).long().squeeze()
        print("Dataset ready.")
            
    def __len__(self):
        return self.X.size(0)
    
    def __getitem__(self, idx):
        # Already normalized and tensorized
        # Add channel dim: (1, 1024, 3)
        return self.X[idx].unsqueeze(0), self.Y[idx]

# --- Model (Exact copy of CNNLSTM's CNN part) ---
class SimpleCNN(nn.Module):
    def __init__(self, num_classes=8):
        super(SimpleCNN, self).__init__()
        
        # Exact architecture from train_frequency_predictor_colab.py
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
        
        # Calculate output size
        # Input: (1, 1024, 3)
        # L1 Pool (4,1): 1024->256
        # L2 Pool (4,1): 256->64
        # Final shape: (32, 64, 3) -> 32*64*3 features
        self.flat_size = 32 * 64 * 3 
        
        self.fc = nn.Linear(self.flat_size, num_classes)

    def forward(self, x):
        # x: (Batch, 1, 1024, 3)
        out = self.cnn(x)
        out = out.view(out.size(0), -1)
        out = self.fc(out)
        return out

def train():
    print("--- Starting CNN Sanity Check ---")
    
    # Path to RAW dataset (not sequences)
    dataset_path = os.path.join('data', 'synthetic', 'classification_dataset_stft_random_deg10_python.h5')
    
    if not os.path.exists(dataset_path):
        print(f"Error: {dataset_path} not found.")
        return

    # Create Dataset & Split
    full_dataset = SingleFrameDataset(dataset_path)
    train_size = int(0.8 * len(full_dataset))
    val_size = len(full_dataset) - train_size
    train_dataset, val_dataset = torch.utils.data.random_split(full_dataset, [train_size, val_size])
    
    train_loader = DataLoader(train_dataset, batch_size=BATCH_SIZE, shuffle=True, num_workers=0)
    val_loader = DataLoader(val_dataset, batch_size=BATCH_SIZE, shuffle=False, num_workers=0)
    
    print(f"Train samples: {len(train_dataset)}, Val samples: {len(val_dataset)}")
    
    # Init Model
    model = SimpleCNN(num_classes=8).to(DEVICE)
    criterion = nn.CrossEntropyLoss()
    optimizer = optim.Adam(model.parameters(), lr=LEARNING_RATE)
    
    # Training Loop
    for epoch in range(EPOCHS):
        model.train()
        train_loss = 0
        correct = 0
        total = 0
        
        for inputs, labels in train_loader:
            inputs, labels = inputs.to(DEVICE), labels.to(DEVICE)
            
            optimizer.zero_grad()
            outputs = model(inputs)
            loss = criterion(outputs, labels)
            loss.backward()
            optimizer.step()
            
            train_loss += loss.item()
            _, predicted = outputs.max(1)
            total += labels.size(0)
            correct += predicted.eq(labels).sum().item()
            
        train_acc = 100. * correct / total
        avg_train_loss = train_loss / len(train_loader)
        
        # Validation
        model.eval()
        val_loss = 0
        correct = 0
        total = 0
        with torch.no_grad():
            for inputs, labels in val_loader:
                inputs, labels = inputs.to(DEVICE), labels.to(DEVICE)
                outputs = model(inputs)
                loss = criterion(outputs, labels)
                
                val_loss += loss.item()
                _, predicted = outputs.max(1)
                total += labels.size(0)
                correct += predicted.eq(labels).sum().item()
        
        val_acc = 100. * correct / total
        avg_val_loss = val_loss / len(val_loader)
        
        print(f"Epoch {epoch+1}/{EPOCHS} | Loss: {avg_train_loss:.4f} | Acc: {train_acc:.2f}% | Val Loss: {avg_val_loss:.4f} | Val Acc: {val_acc:.2f}%")
    
    # Save the trained weights
    if not os.path.exists('models'):
        os.makedirs('models')
    weight_path = os.path.join('models', 'cnn_weights.pth')
    torch.save(model.cnn.state_dict(), weight_path)
    print(f"Saved CNN weights to {weight_path}")

if __name__ == "__main__":
    # Fix for python 2 long
    if 'long' not in locals():
        long = int
    train()
