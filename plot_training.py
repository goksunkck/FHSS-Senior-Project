import matplotlib.pyplot as plt
import re

# Save your training output to a file first, or read from the cell output
# Option 1: If you saved output to a file
# with open('training_log.txt', 'r') as f:
#     lines = f.readlines()

# Option 2: Paste the log into this variable
log = open('training_log.txt').readlines()

epochs, train_loss, val_loss, bit_acc, state_acc = [], [], [], [], []

for line in log:
    m = re.search(r'Epoch (\d+)/\d+ \| Train Loss: ([\d.]+) \| Val Loss: ([\d.]+) \| Val Bit Acc: ([\d.]+)% \| Exact 20-Bit State Acc: ([\d.]+)%', line)
    if m:
        epochs.append(int(m.group(1)))
        train_loss.append(float(m.group(2)))
        val_loss.append(float(m.group(3)))
        bit_acc.append(float(m.group(4)))
        state_acc.append(float(m.group(5)))

fig, axes = plt.subplots(1, 3, figsize=(20, 5))

axes[0].plot(epochs, train_loss, label='Train Loss (CoT)', color='blue')
axes[0].plot(epochs, val_loss, label='Val Loss (Autoregressive)', color='orange')
axes[0].set_xlabel('Epoch'); axes[0].set_ylabel('MSE Loss')
axes[0].set_title('Loss Curves'); axes[0].legend(); axes[0].grid(True, alpha=0.3)

axes[1].plot(epochs, bit_acc, label='Val Bit Accuracy', color='green')
axes[1].set_xlabel('Epoch'); axes[1].set_ylabel('Accuracy (%)')
axes[1].set_title('Per-Bit Accuracy'); axes[1].legend(); axes[1].grid(True, alpha=0.3)

axes[2].plot(epochs, state_acc, label='Exact 20-Bit State Acc', color='red')
axes[2].set_xlabel('Epoch'); axes[2].set_ylabel('Accuracy (%)')
axes[2].set_title('Exact 20-Bit State Accuracy'); axes[2].legend(); axes[2].grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('models/state_transformer_training.png', dpi=150, bbox_inches='tight')
plt.show()
