import argparse
import math
import os
import sys

import h5py
import matplotlib.pyplot as plt
import numpy as np
import torch
import torch.nn as nn
import torch.optim as optim
import torch.nn.functional as F
from torch.utils.data import DataLoader, Dataset

sys.path.append(os.path.join(os.getcwd(), "src", "python"))
from dsp_utils import add_awgn, compute_stft_matrix, fsk_modulate


DEVICE = torch.device("cuda" if torch.cuda.is_available() else "cpu")


class MSequenceLFSR:
    def __init__(self, initial_state, taps=(10, 3)):
        self.state = np.array(initial_state, dtype=np.uint8)
        self.taps = taps

    def step(self):
        feedback = 0
        for tap in self.taps:
            feedback ^= int(self.state[tap - 1])
        output = int(self.state[-1])
        self.state = np.roll(self.state, 1)
        self.state[0] = feedback
        return output

    def generate(self, n_bits):
        return np.array([self.step() for _ in range(n_bits)], dtype=np.uint8)

    def state_bits(self):
        return self.state.astype(np.float32)


def seed_to_state(seed_int):
    return [int(bit) for bit in format(int(seed_int), "010b")]


def state_to_hop(bits):
    arr = np.asarray(bits, dtype=np.uint8).reshape(-1)
    lfsr = MSequenceLFSR(arr)
    hop_bits = np.array([lfsr.step(), lfsr.step(), lfsr.step()], dtype=np.uint8)
    return int(hop_bits.dot(np.array([4, 2, 1], dtype=np.uint8)))


def generate_m_sequence_dataset(output_path, total_signals, seed):
    print(f"Generating compact m-sequence dataset: {total_signals} signals -> {output_path}")
    os.makedirs(os.path.dirname(output_path), exist_ok=True)

    rng = np.random.default_rng(seed)
    snr_levels_db = np.arange(-10, 12, 2)
    fs = 10e6
    num_bits = 512
    symbol_rate = fs * 0.01
    freq_separation = symbol_rate
    samples_per_symbol = int(round(fs / symbol_rate))
    k = 3
    num_hops = 256
    num_channels = 2**k
    spacing = 2 * freq_separation
    base_freq = 2e6
    hopset = np.arange(num_channels) * spacing + base_freq
    sim_params = {
        "fs": fs,
        "M": 2,
        "freqSeparation": freq_separation,
        "samplesPerSymbol": samples_per_symbol,
    }
    stft_win = 128
    stft_over = 96
    stft_nfft = 1024
    samples_per_hop = int(np.floor((num_bits * samples_per_symbol) / num_hops))

    all_seeds = np.arange(1, 1024)
    rng.shuffle(all_seeds)
    train_seed_end = int(0.70 * len(all_seeds))
    val_seed_end = train_seed_end + int(0.15 * len(all_seeds))
    seed_pools = [
        all_seeds[:train_seed_end],
        all_seeds[train_seed_end:val_seed_end],
        all_seeds[val_seed_end:],
    ]
    split_counts = [
        int(0.70 * total_signals),
        int(0.15 * total_signals),
        total_signals - int(0.70 * total_signals) - int(0.15 * total_signals),
    ]

    total_examples = total_signals * num_hops
    global_min = float("inf")
    global_max = float("-inf")
    write_idx = 0

    with h5py.File(output_path, "w") as f:
        dset_x = f.create_dataset(
            "X_data",
            (total_examples, 1024, 3),
            dtype="float32",
            chunks=(min(64 * num_hops, total_examples), 1024, 3),
        )
        dset_y = f.create_dataset("Y_data", (total_examples, 1), dtype="uint8")
        dset_set = f.create_dataset("Set_ID", (total_examples, 1), dtype="uint8")
        dset_seed = f.create_dataset("Seed_Val", (total_examples, 1), dtype="int32")
        dset_snr = f.create_dataset("SNR", (total_examples, 1), dtype="float32")
        f.attrs["numHopsPerSignal"] = num_hops
        f.attrs["polynomial_degree"] = 10
        f.attrs["sequence_type"] = "m_sequence"

        for set_id, (pool, count) in enumerate(zip(seed_pools, split_counts)):
            for i in range(count):
                seed_int = int(rng.choice(pool))
                snr_db = float(rng.choice(snr_levels_db))
                lfsr = MSequenceLFSR(seed_to_state(seed_int))
                seq_bits = lfsr.generate(num_hops * k)
                hop_bits = seq_bits.reshape(k, -1, order="F").T
                hop_indices = hop_bits.dot(np.array([4, 2, 1], dtype=np.uint8)).astype(np.uint8)
                hop_freqs = hopset[hop_indices]

                message_bits = rng.integers(0, 2, num_bits, dtype=np.uint8)
                modulated_data = fsk_modulate(message_bits, sim_params)
                t_vec = np.arange(samples_per_hop) / fs
                x_batch = np.zeros((num_hops, 1024, 3), dtype=np.float32)

                for hop in range(num_hops):
                    start = hop * samples_per_hop
                    end = (hop + 1) * samples_per_hop
                    carrier = np.exp(1j * 2 * np.pi * hop_freqs[hop] * t_vec)
                    rx_signal = add_awgn(modulated_data[start:end] * carrier, snr_db)
                    stft_log = compute_stft_matrix(rx_signal, fs, stft_nfft, stft_win, stft_over)
                    x_batch[hop] = stft_log

                global_min = min(global_min, float(np.min(x_batch)))
                global_max = max(global_max, float(np.max(x_batch)))

                end_idx = write_idx + num_hops
                dset_x[write_idx:end_idx] = x_batch
                dset_y[write_idx:end_idx, 0] = hop_indices
                dset_set[write_idx:end_idx, 0] = set_id
                dset_seed[write_idx:end_idx, 0] = seed_int
                dset_snr[write_idx:end_idx, 0] = snr_db
                write_idx = end_idx

                if (i + 1) % 10 == 0 or i + 1 == count:
                    print(f"  split {set_id}: {i + 1}/{count} signals")

        f.create_dataset("global_min", data=np.array([global_min], dtype=np.float32))
        f.create_dataset("global_max", data=np.array([global_max], dtype=np.float32))
    print(f"Dataset complete. Global min={global_min:.2f}, max={global_max:.2f}")


def build_sequence_starts(h5_path, lookback_window):
    with h5py.File(h5_path, "r") as f:
        set_ids = f["Set_ID"][:, 0]
        num_hops = int(f.attrs.get("numHopsPerSignal", 256))

    starts_by_split = {}
    for set_id in (0, 1, 2):
        signal_offsets = np.flatnonzero(set_ids[::num_hops] == set_id) * num_hops
        starts = [
            offset + local
            for offset in signal_offsets
            for local in range(num_hops - lookback_window)
        ]
        starts_by_split[set_id] = np.asarray(starts, dtype=np.int64)
    return starts_by_split


class MSequenceStateDataset(Dataset):
    def __init__(self, h5_path, sequence_starts, lookback_window, load_ram=True, norm_min=None, norm_max=None):
        self.h5_path = h5_path
        self.sequence_starts = np.asarray(sequence_starts, dtype=np.int64)
        self.lookback_window = lookback_window
        self.load_ram = load_ram
        self.file_handle = None

        with h5py.File(h5_path, "r") as f:
            self.seed_vals = f["Seed_Val"][:, 0].astype(np.int32)
            self.next_hops = f["Y_data"][:, 0].astype(np.int64)
            self.num_hops = int(f.attrs.get("numHopsPerSignal", 256))
            file_min = float(f["global_min"][0]) if "global_min" in f else -120.0
            file_max = float(f["global_max"][0]) if "global_max" in f else 0.0
            self.global_min = file_min if norm_min is None else float(norm_min)
            self.global_max = file_max if norm_max is None else float(norm_max)
            if load_ram:
                print(f"Loading {h5_path} into RAM...")
                x_data = f["X_data"][:]
                x_data = (x_data - self.global_min) / (self.global_max - self.global_min + 1e-6)
                self.x_data = torch.from_numpy(x_data.astype(np.float32, copy=False)).unsqueeze(1)
            else:
                self.x_data = None

        print("Precomputing 10-bit m-sequence states...")
        self.state_cache = {}
        for seed_int in range(1, 1024):
            lfsr = MSequenceLFSR(seed_to_state(seed_int))
            states = []
            for _ in range(self.num_hops):
                bipolar = lfsr.state_bits() * 2.0 - 1.0
                states.append(torch.tensor(bipolar, dtype=torch.float32))
                for _ in range(3):
                    lfsr.step()
            self.state_cache[seed_int] = states
        self.state_cache[0] = [torch.zeros(10, dtype=torch.float32) for _ in range(self.num_hops)]

    def __len__(self):
        return len(self.sequence_starts)

    def _h5(self):
        if self.file_handle is None:
            self.file_handle = h5py.File(self.h5_path, "r")
        return self.file_handle

    def __getitem__(self, idx):
        start = int(self.sequence_starts[idx])
        end = start + self.lookback_window
        target = end

        if self.load_ram:
            seq = self.x_data[start:end]
        else:
            raw = self._h5()["X_data"][start:end]
            raw = (raw - self.global_min) / (self.global_max - self.global_min + 1e-6)
            seq = torch.from_numpy(raw.astype(np.float32, copy=False)).unsqueeze(1)

        local_target_hop = target % self.num_hops
        seed_int = int(self.seed_vals[target])
        state = self.state_cache[seed_int][local_target_hop]
        next_hop = torch.tensor(self.next_hops[target], dtype=torch.long)
        return seq, state, next_hop


class FrameDataset(Dataset):
    def __init__(self, h5_path, indices, norm_min=None, norm_max=None):
        self.h5_path = h5_path
        self.indices = np.asarray(indices, dtype=np.int64)
        self.file_handle = None
        with h5py.File(h5_path, "r") as f:
            file_min = float(f["global_min"][0]) if "global_min" in f else -120.0
            file_max = float(f["global_max"][0]) if "global_max" in f else 0.0
            self.global_min = file_min if norm_min is None else float(norm_min)
            self.global_max = file_max if norm_max is None else float(norm_max)

    def _h5(self):
        if self.file_handle is None:
            self.file_handle = h5py.File(self.h5_path, "r")
        return self.file_handle

    def __len__(self):
        return len(self.indices)

    def __getitem__(self, idx):
        source_idx = int(self.indices[idx])
        f = self._h5()
        x = f["X_data"][source_idx]
        x = (x - self.global_min) / (self.global_max - self.global_min + 1e-6)
        y = int(f["Y_data"][source_idx, 0])
        return torch.from_numpy(x.astype(np.float32, copy=False)).unsqueeze(0), torch.tensor(y, dtype=torch.long)


class PositionalEncoding(nn.Module):
    def __init__(self, d_model, max_len=5000):
        super().__init__()
        pe = torch.zeros(max_len, d_model)
        position = torch.arange(0, max_len, dtype=torch.float).unsqueeze(1)
        div_term = torch.exp(torch.arange(0, d_model, 2).float() * (-math.log(10000.0) / d_model))
        pe[:, 0::2] = torch.sin(position * div_term)
        pe[:, 1::2] = torch.cos(position * div_term)
        pe = pe.unsqueeze(0).transpose(0, 1)
        self.register_buffer("pe", pe)

    def forward(self, x):
        return x + self.pe[: x.size(0), :]


class CosineActivation(nn.Module):
    def forward(self, x):
        return torch.cos(math.pi * x)


class StateEncoderLayer(nn.TransformerEncoderLayer):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.activation = CosineActivation()


class MSequenceStateTransformer(nn.Module):
    def __init__(self, state_dim=10, num_classes=8, d_model=64, nhead=4, num_layers=2):
        super().__init__()
        self.state_dim = state_dim
        self.num_classes = num_classes
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
            nn.ReLU(),
        )
        self.flat_size = 32 * 64 * 3
        self.cnn_fc = nn.Linear(self.flat_size, num_classes)
        self.prob_proj = nn.Linear(num_classes, d_model)
        self.bit_proj = nn.Linear(1, d_model)
        self.pos_encoder = PositionalEncoding(d_model)
        encoder_layer = StateEncoderLayer(
            d_model=d_model,
            nhead=nhead,
            dim_feedforward=128,
            dropout=0.0,
        )
        self.transformer_encoder = nn.TransformerEncoder(encoder_layer, num_layers=num_layers)
        self.fc = nn.Linear(d_model, 1)

    def generate_square_subsequent_mask(self, size):
        mask = (torch.triu(torch.ones(size, size)) == 1).transpose(0, 1)
        return mask.float().masked_fill(mask == 0, float("-inf")).masked_fill(mask == 1, float(0.0))

    def visual_tokens(self, x):
        batch, seq_len, channels, height, width = x.size()
        c_in = x.view(batch * seq_len, channels, height, width)
        with torch.no_grad():
            c_out = self.cnn(c_in)
            logits = self.cnn_fc(c_out.view(c_out.size(0), -1))
            preds = torch.argmax(logits, dim=-1)
            one_hot = F.one_hot(preds, num_classes=self.num_classes).float()
        return self.prob_proj(one_hot).view(batch, seq_len, -1)

    def forward(self, x, y_state=None):
        batch, seq_len = x.size(0), x.size(1)
        vis_emb = self.visual_tokens(x)

        if y_state is not None:
            teacher_bits = y_state[:, :-1].unsqueeze(-1)
            state_emb = self.bit_proj(teacher_bits)
            src = torch.cat([vis_emb, state_emb], dim=1)
            src = src.permute(1, 0, 2)
            src = self.pos_encoder(src)
            mask = self.generate_square_subsequent_mask(src.size(0)).to(src.device)
            out = self.transformer_encoder(src, mask=mask)
            state_preds = self.fc(out[seq_len - 1 : seq_len - 1 + self.state_dim, :, :])
            return state_preds.squeeze(-1).permute(1, 0)

        curr_src = vis_emb
        generated_bits = []
        for step in range(self.state_dim):
            src = curr_src.permute(1, 0, 2)
            src = self.pos_encoder(src)
            mask = self.generate_square_subsequent_mask(src.size(0)).to(src.device)
            out = self.transformer_encoder(src, mask=mask)
            next_bit_logit = self.fc(out[-1, :, :])
            generated_bits.append(next_bit_logit)
            if step < self.state_dim - 1:
                next_bit = (next_bit_logit > 0).float() * 2.0 - 1.0
                curr_src = torch.cat([curr_src, self.bit_proj(next_bit).unsqueeze(1)], dim=1)
        return torch.cat(generated_bits, dim=-1)


def load_cnn_weights(model, weight_path):
    if not weight_path or not os.path.exists(weight_path):
        print("CNN weights not found; the visual front-end will remain randomly initialized.")
        return
    print(f"Loading CNN weights from {weight_path}")
    state = torch.load(weight_path, map_location=DEVICE)
    cnn_state = {}
    fc_state = {}
    for key, value in state.items():
        if key.startswith("cnn."):
            cnn_state[key[4:]] = value
        elif key.startswith("fc."):
            fc_state[key[3:]] = value
        elif key[0].isdigit():
            cnn_state[key] = value
        elif key in ("weight", "bias"):
            fc_state[key] = value
    model.cnn.load_state_dict(cnn_state, strict=True)
    if fc_state.get("weight") is not None and fc_state["weight"].shape == model.cnn_fc.weight.shape:
        model.cnn_fc.weight.data.copy_(fc_state["weight"])
        model.cnn_fc.bias.data.copy_(fc_state["bias"])
    print("CNN weights loaded.")


def pretrain_cnn(model, h5_path, args, norm_min=None, norm_max=None):
    print(f"Pretraining CNN front-end on m-sequence frames for {args.pretrain_cnn_epochs} epoch(s)...")
    with h5py.File(h5_path, "r") as f:
        set_ids = f["Set_ID"][:, 0]
    train_indices = np.flatnonzero(set_ids == 0)
    val_indices = np.flatnonzero(set_ids == 1)
    rng = np.random.default_rng(args.seed)
    rng.shuffle(train_indices)
    rng.shuffle(val_indices)
    if args.max_cnn_train:
        train_indices = train_indices[: args.max_cnn_train]
    if args.max_cnn_val:
        val_indices = val_indices[: args.max_cnn_val]

    train_ds = FrameDataset(h5_path, train_indices, norm_min=norm_min, norm_max=norm_max)
    val_ds = FrameDataset(h5_path, val_indices, norm_min=norm_min, norm_max=norm_max)
    train_loader = DataLoader(train_ds, batch_size=args.batch_size, shuffle=True, num_workers=0)
    val_loader = DataLoader(val_ds, batch_size=args.batch_size, shuffle=False, num_workers=0)

    for param in model.cnn.parameters():
        param.requires_grad = True
    for param in model.cnn_fc.parameters():
        param.requires_grad = True

    optimizer = optim.Adam(list(model.cnn.parameters()) + list(model.cnn_fc.parameters()), lr=args.cnn_lr)
    criterion = nn.CrossEntropyLoss()
    for epoch in range(args.pretrain_cnn_epochs):
        model.train()
        train_loss = 0.0
        train_correct = 0
        train_total = 0
        for x, y in train_loader:
            x = x.to(DEVICE)
            y = y.to(DEVICE)
            optimizer.zero_grad()
            logits = model.cnn_fc(model.cnn(x).view(x.size(0), -1))
            loss = criterion(logits, y)
            loss.backward()
            optimizer.step()
            train_loss += loss.item()
            train_correct += (logits.argmax(1) == y).sum().item()
            train_total += y.size(0)

        model.eval()
        val_correct = 0
        val_total = 0
        val_loss = 0.0
        with torch.no_grad():
            for x, y in val_loader:
                x = x.to(DEVICE)
                y = y.to(DEVICE)
                logits = model.cnn_fc(model.cnn(x).view(x.size(0), -1))
                val_loss += criterion(logits, y).item()
                val_correct += (logits.argmax(1) == y).sum().item()
                val_total += y.size(0)

        print(
            f"CNN Epoch {epoch + 1}/{args.pretrain_cnn_epochs} | "
            f"Loss {train_loss / max(1, len(train_loader)):.4f} | "
            f"Train Acc {100.0 * train_correct / max(1, train_total):.1f}% | "
            f"Val Loss {val_loss / max(1, len(val_loader)):.4f} | "
            f"Val Acc {100.0 * val_correct / max(1, val_total):.1f}%"
        )

    if args.mseq_cnn_checkpoint:
        os.makedirs(os.path.dirname(args.mseq_cnn_checkpoint), exist_ok=True)
        torch.save(
            {
                **{f"cnn.{key}": value.detach().cpu() for key, value in model.cnn.state_dict().items()},
                **{f"fc.{key}": value.detach().cpu() for key, value in model.cnn_fc.state_dict().items()},
            },
            args.mseq_cnn_checkpoint,
        )
        print(f"saved m-sequence CNN weights: {args.mseq_cnn_checkpoint}")


def predicted_hops_from_states(pred_bits):
    bits = (pred_bits.detach().cpu().numpy() > 0).astype(np.uint8)
    return np.array([state_to_hop(row) for row in bits], dtype=np.int64)


def train(args):
    if args.plot_only:
        checkpoint = torch.load(args.checkpoint, map_location="cpu")
        if not isinstance(checkpoint, dict) or "history" not in checkpoint:
            raise ValueError(f"No training history found in {args.checkpoint}")
        save_training_plot(checkpoint["history"], args.plot)
        return

    if not os.path.exists(args.dataset):
        if not args.generate_if_missing:
            raise FileNotFoundError(
                f"{args.dataset} is missing. Re-run with --generate-if-missing or provide --dataset."
            )
        generate_m_sequence_dataset(args.dataset, args.signals, args.seed)

    starts = build_sequence_starts(args.dataset, args.lookback)
    rng = np.random.default_rng(args.seed)
    train_starts = starts[0]
    val_starts = starts[1]
    rng.shuffle(train_starts)
    rng.shuffle(val_starts)
    if args.max_train:
        train_starts = train_starts[: args.max_train]
    if args.max_val:
        val_starts = val_starts[: args.max_val]

    norm_min = args.norm_min
    norm_max = args.norm_max
    print(f"Train windows: {len(train_starts)}, validation windows: {len(val_starts)}")

    model = MSequenceStateTransformer(state_dim=10).to(DEVICE)
    load_cnn_weights(model, args.cnn_weights)

    history = {"train_loss": [], "val_loss": [], "val_bit_acc": [], "val_state_acc": [], "val_hop_acc": []}
    if args.resume and os.path.exists(args.checkpoint):
        checkpoint = torch.load(args.checkpoint, map_location=DEVICE)
        state_dict = checkpoint["model_state_dict"] if isinstance(checkpoint, dict) and "model_state_dict" in checkpoint else checkpoint
        model.load_state_dict(state_dict, strict=True)
        if isinstance(checkpoint, dict) and isinstance(checkpoint.get("history"), dict):
            history.update(checkpoint["history"])
        print(f"Resumed state transformer from {args.checkpoint}")

    if args.pretrain_cnn_epochs > 0:
        pretrain_cnn(model, args.dataset, args, norm_min=norm_min, norm_max=norm_max)

    train_ds = MSequenceStateDataset(
        args.dataset,
        train_starts,
        args.lookback,
        load_ram=not args.no_load_ram,
        norm_min=norm_min,
        norm_max=norm_max,
    )
    val_ds = MSequenceStateDataset(
        args.dataset,
        val_starts,
        args.lookback,
        load_ram=not args.no_load_ram,
        norm_min=norm_min,
        norm_max=norm_max,
    )
    train_loader = DataLoader(train_ds, batch_size=args.batch_size, shuffle=True, num_workers=0)
    val_loader = DataLoader(val_ds, batch_size=args.batch_size, shuffle=False, num_workers=0)

    for param in model.cnn.parameters():
        param.requires_grad = False
    for param in model.cnn_fc.parameters():
        param.requires_grad = False

    optimizer = optim.AdamW(filter(lambda p: p.requires_grad, model.parameters()), lr=args.lr)
    scheduler = optim.lr_scheduler.CosineAnnealingLR(optimizer, T_max=args.epochs)
    criterion = nn.MSELoss()
    best_val_loss = min(history["val_loss"]) if history["val_loss"] else float("inf")
    best_val_hop_acc = max(history["val_hop_acc"]) if history["val_hop_acc"] else float("-inf")
    start_epoch = len(history["train_loss"])

    for epoch in range(args.epochs):
        model.train()
        train_loss = 0.0
        train_state_correct = 0
        train_total = 0
        for x, y_state, _ in train_loader:
            x = x.to(DEVICE)
            y_state = y_state.to(DEVICE)
            optimizer.zero_grad()
            out = model(x, y_state=y_state)
            loss = criterion(out, y_state)
            loss.backward()
            optimizer.step()
            train_loss += loss.item()
            pred_bits = (out > 0).float() * 2.0 - 1.0
            train_state_correct += (pred_bits == y_state).all(dim=1).sum().item()
            train_total += y_state.size(0)

        model.eval()
        val_loss = 0.0
        val_bit_correct = 0
        val_state_correct = 0
        val_hop_correct = 0
        val_total = 0
        with torch.no_grad():
            for x, y_state, y_hop in val_loader:
                x = x.to(DEVICE)
                y_state = y_state.to(DEVICE)
                out = model(x)
                val_loss += criterion(out, y_state).item()
                pred_bits = (out > 0).float() * 2.0 - 1.0
                val_bit_correct += (pred_bits == y_state).sum().item()
                val_state_correct += (pred_bits == y_state).all(dim=1).sum().item()
                pred_hops = predicted_hops_from_states(out)
                val_hop_correct += int(np.sum(pred_hops == y_hop.numpy()))
                val_total += y_state.size(0)

        avg_train_loss = train_loss / max(1, len(train_loader))
        avg_val_loss = val_loss / max(1, len(val_loader))
        train_state_acc = 100.0 * train_state_correct / max(1, train_total)
        val_bit_acc = 100.0 * val_bit_correct / max(1, val_total * 10)
        val_state_acc = 100.0 * val_state_correct / max(1, val_total)
        val_hop_acc = 100.0 * val_hop_correct / max(1, val_total)
        scheduler.step()

        history["train_loss"].append(avg_train_loss)
        history["val_loss"].append(avg_val_loss)
        history["val_bit_acc"].append(val_bit_acc)
        history["val_state_acc"].append(val_state_acc)
        history["val_hop_acc"].append(val_hop_acc)

        print(
            f"Epoch {start_epoch + epoch + 1} (+{epoch + 1}/{args.epochs}) | "
            f"Train Loss {avg_train_loss:.4f} | "
            f"Train State {train_state_acc:.1f}% | Val Loss {avg_val_loss:.4f} | "
            f"Val Bit {val_bit_acc:.1f}% | Val State {val_state_acc:.1f}% | Val Hop {val_hop_acc:.1f}%"
        )

        should_save = avg_val_loss < best_val_loss or val_hop_acc > best_val_hop_acc
        if should_save:
            best_val_loss = avg_val_loss
            best_val_hop_acc = max(best_val_hop_acc, val_hop_acc)
            os.makedirs(os.path.dirname(args.checkpoint), exist_ok=True)
            torch.save(
                {
                    "model_state_dict": model.state_dict(),
                    "lookback": args.lookback,
                    "state_dim": 10,
                    "dataset": args.dataset,
                    "history": history,
                },
                args.checkpoint,
            )
            print(f"  saved checkpoint: {args.checkpoint}")

        if args.target_hop_acc and val_hop_acc >= args.target_hop_acc:
            print(f"Reached target hop accuracy: {val_hop_acc:.3f}% >= {args.target_hop_acc:.3f}%")
            break
        if args.target_state_acc and val_state_acc >= args.target_state_acc:
            print(f"Reached target state accuracy: {val_state_acc:.3f}% >= {args.target_state_acc:.3f}%")
            break

    if args.plot:
        save_training_plot(history, args.plot)


def save_training_plot(history, out_path):
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    epochs = np.arange(1, len(history["train_loss"]) + 1)
    fig, axes = plt.subplots(1, 3, figsize=(18, 4.5), dpi=160)

    axes[0].plot(epochs, history["train_loss"], color="blue", linewidth=1.8, label="Train Loss")
    axes[0].plot(epochs, history["val_loss"], color="orange", linewidth=1.8, label="Val Loss")
    axes[0].set_title("Loss Curves")
    axes[0].set_xlabel("Epoch")
    axes[0].set_ylabel("MSE Loss")
    axes[0].grid(True, alpha=0.3)
    axes[0].legend(loc="upper right")

    axes[1].plot(epochs, history["val_bit_acc"], color="green", linewidth=1.8, label="Val Bit Accuracy")
    axes[1].set_title("Per-Bit Accuracy")
    axes[1].set_xlabel("Epoch")
    axes[1].set_ylabel("Accuracy (%)")
    axes[1].set_ylim(0, 102)
    axes[1].grid(True, alpha=0.3)
    axes[1].legend(loc="lower right")

    axes[2].plot(epochs, history["val_state_acc"], color="red", linewidth=1.8, label="Exact 10-Bit State Acc")
    axes[2].set_title("Exact 10-Bit State Accuracy")
    axes[2].set_xlabel("Epoch")
    axes[2].set_ylabel("Accuracy (%)")
    axes[2].set_ylim(0, 102)
    axes[2].grid(True, alpha=0.3)
    axes[2].legend(loc="lower right")

    fig.tight_layout()
    fig.savefig(out_path, bbox_inches="tight")
    plt.close(fig)
    print(f"saved plot: {out_path}")


def parse_args():
    parser = argparse.ArgumentParser(description="Train the state transformer on m-sequence internal states.")
    parser.add_argument("--dataset", default=os.path.join("data", "synthetic", "m_sequence_state_dataset.h5"))
    parser.add_argument("--checkpoint", default=os.path.join("models", "m_sequence_state_transformer_best.pth"))
    parser.add_argument("--cnn-weights", default=os.path.join("models", "cnn_weights.pth"))
    parser.add_argument("--mseq-cnn-checkpoint", default=os.path.join("models", "m_sequence_cnn_weights.pth"))
    parser.add_argument("--plot", default=os.path.join("results", "m_sequence_state_transformer_training.png"))
    parser.add_argument("--generate-if-missing", action="store_true")
    parser.add_argument("--signals", type=int, default=120)
    parser.add_argument("--epochs", type=int, default=1)
    parser.add_argument("--pretrain-cnn-epochs", type=int, default=0)
    parser.add_argument("--resume", action="store_true")
    parser.add_argument("--plot-only", action="store_true")
    parser.add_argument("--target-hop-acc", type=float, default=0.0)
    parser.add_argument("--target-state-acc", type=float, default=0.0)
    parser.add_argument("--batch-size", type=int, default=128)
    parser.add_argument("--lr", type=float, default=1e-3)
    parser.add_argument("--cnn-lr", type=float, default=1e-3)
    parser.add_argument("--lookback", type=int, default=12)
    parser.add_argument("--max-train", type=int, default=0)
    parser.add_argument("--max-val", type=int, default=0)
    parser.add_argument("--max-cnn-train", type=int, default=0)
    parser.add_argument("--max-cnn-val", type=int, default=0)
    parser.add_argument("--norm-min", type=float, default=None)
    parser.add_argument("--norm-max", type=float, default=None)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--no-load-ram", action="store_true")
    return parser.parse_args()


if __name__ == "__main__":
    train(parse_args())
