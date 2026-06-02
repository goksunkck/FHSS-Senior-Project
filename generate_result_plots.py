import argparse
import json
import os
import sys

import h5py
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import torch

sys.path.append(os.path.join(os.getcwd(), "src", "python"))

from berlekamp_massey import berlekamp_massey, reconstruct_lfsr_state
from gold_code import LFSR
from train_m_sequence_state_transformer import (
    MSequenceStateTransformer,
    predicted_hops_from_states,
)
from train_state_transformer import StateTransformer


DEVICE = torch.device("cuda" if torch.cuda.is_available() else "cpu")
LOOKBACK = 12
NUM_HOPS = 256
BITS_PER_HOP = 3


def require_file(path, description):
    if not os.path.exists(path):
        raise FileNotFoundError(f"{description} not found: {path}")


def ensure_dir(path):
    os.makedirs(path, exist_ok=True)


def load_checkpoint(path):
    checkpoint = torch.load(path, map_location=DEVICE)
    if isinstance(checkpoint, dict) and "model_state_dict" in checkpoint:
        return checkpoint["model_state_dict"], checkpoint
    return checkpoint, {}


def load_history(path):
    if not path or not os.path.exists(path):
        return None
    if path.lower().endswith(".json"):
        with open(path, "r", encoding="utf-8") as fp:
            return json.load(fp)
    checkpoint = torch.load(path, map_location="cpu")
    if isinstance(checkpoint, dict):
        return checkpoint.get("history")
    return None


def load_signal(h5_path, signal_index):
    with h5py.File(h5_path, "r") as f:
        base = signal_index * NUM_HOPS
        x_raw = f["X_data"][base : base + NUM_HOPS]
        y = f["Y_data"][base : base + NUM_HOPS, 0].astype(int)
        snr = float(f["SNR"][base, 0]) if "SNR" in f else float("nan")
        seed = int(f["Seed_Val"][base, 0]) if "Seed_Val" in f else -1

        if "global_min" in f and "global_max" in f:
            x_min = float(f["global_min"][0])
            x_max = float(f["global_max"][0])
        else:
            x_min = float(np.min(x_raw))
            x_max = float(np.max(x_raw))

    x = (x_raw - x_min) / (x_max - x_min + 1e-6)
    return x.astype(np.float32), y, snr, seed


def state_to_next_hop(state_bits_20):
    binary_bits = [1 if bit > 0 else 0 for bit in state_bits_20]

    class StateLFSR:
        def __init__(self, taps, state_bits):
            self.degree = 10
            self.taps = taps
            self.state = sum(int(bit) << i for i, bit in enumerate(state_bits))

        def step(self):
            out_bit = (self.state >> (self.degree - 1)) & 1
            feedback = out_bit
            for tap in self.taps:
                feedback ^= (self.state >> (tap - 1)) & 1
            self.state = ((self.state << 1) & ((1 << self.degree) - 1)) | feedback
            return out_bit

    lfsr1 = StateLFSR([3], binary_bits[:10])
    lfsr2 = StateLFSR([7], binary_bits[10:])
    chips = [lfsr1.step() ^ lfsr2.step() for _ in range(BITS_PER_HOP)]
    return chips[0] * 4 + chips[1] * 2 + chips[2]


def predict_gold_signal(model, x_signal, batch_size=64):
    x_tensor = torch.from_numpy(x_signal).float().unsqueeze(1).to(DEVICE)
    windows = torch.stack(
        [x_tensor[i : i + LOOKBACK] for i in range(NUM_HOPS - LOOKBACK)],
        dim=0,
    )

    pred_hops = []
    with torch.no_grad():
        cnn_logits = model.cnn_fc(model.cnn(x_tensor).view(NUM_HOPS, -1))
        cnn_preds = torch.argmax(cnn_logits, dim=1).cpu().numpy()

        for start in range(0, len(windows), batch_size):
            state_out = model(windows[start : start + batch_size])
            pred_hops.extend(state_to_next_hop(row) for row in state_out.cpu().numpy())

    return cnn_preds, np.asarray(pred_hops, dtype=int)


def predict_mseq_signal(model, x_signal, batch_size=64):
    x_tensor = torch.from_numpy(x_signal).float().unsqueeze(1).to(DEVICE)
    windows = torch.stack(
        [x_tensor[i : i + LOOKBACK] for i in range(NUM_HOPS - LOOKBACK)],
        dim=0,
    )

    pred_hops = []
    with torch.no_grad():
        for start in range(0, len(windows), batch_size):
            out = model(windows[start : start + batch_size])
            pred_hops.extend(predicted_hops_from_states(out).tolist())
    return np.asarray(pred_hops, dtype=int)


def bm_predict_from_hops(observed_hops, num_future):
    bit_sequence = []
    for hop in observed_hops:
        bit_sequence.extend([(int(hop) >> (BITS_PER_HOP - 1 - i)) & 1 for i in range(BITS_PER_HOP)])

    polynomial, _ = berlekamp_massey(np.asarray(bit_sequence, dtype=int))
    initial_state = reconstruct_lfsr_state(bit_sequence, polynomial)
    lfsr = LFSR(polynomial, initial_state)
    _ = lfsr.generate(len(bit_sequence))
    future_bits = lfsr.generate(num_future * BITS_PER_HOP)

    predicted = []
    for i in range(0, len(future_bits), BITS_PER_HOP):
        bits = future_bits[i : i + BITS_PER_HOP]
        predicted.append(int(bits[0]) * 4 + int(bits[1]) * 2 + int(bits[2]))
    return np.asarray(predicted, dtype=int)


def plot_spectrogram(x_signal, y_signal, snr, out_path):
    combined = np.concatenate([x_signal[i] for i in range(NUM_HOPS)], axis=1)
    fig, (ax0, ax1) = plt.subplots(
        2,
        1,
        figsize=(14, 7),
        dpi=160,
        gridspec_kw={"height_ratios": [3, 1]},
    )

    im = ax0.imshow(combined, aspect="auto", origin="lower", cmap="viridis", interpolation="nearest")
    ax0.set_title(f"Generated FHSS Spectrogram Example (SNR = {snr:.0f} dB)")
    ax0.set_ylabel("Frequency bin")
    ax0.set_xlabel("STFT time column")
    fig.colorbar(im, ax=ax0, label="Normalized log magnitude")

    ax1.step(np.arange(NUM_HOPS), y_signal, where="mid", linewidth=1.5, color="#1565C0")
    ax1.set_yticks(range(8))
    ax1.set_ylim(-0.5, 7.5)
    ax1.set_xlabel("Hop index")
    ax1.set_ylabel("Channel")
    ax1.grid(True, alpha=0.3)

    fig.tight_layout()
    fig.savefig(out_path, bbox_inches="tight")
    plt.close(fig)


def plot_prediction(actual, predicted, title, out_path, start_hop=LOOKBACK, max_points=120):
    actual = np.asarray(actual)[:max_points]
    predicted = np.asarray(predicted)[:max_points]
    pos = np.arange(start_hop, start_hop + len(actual))
    correct = actual == predicted

    fig, ax = plt.subplots(figsize=(14, 4.5), dpi=160)
    ax.step(pos, actual, where="mid", linewidth=1.8, color="#1565C0", label="Actual hop")
    ax.scatter(pos[correct], predicted[correct], s=22, color="#2E7D32", label="Correct prediction", zorder=4)
    if np.any(~correct):
        ax.scatter(pos[~correct], predicted[~correct], s=45, marker="x", color="#C62828", label="Wrong prediction", zorder=5)
    ax.set_title(title)
    ax.set_xlabel("Hop index")
    ax.set_ylabel("Frequency channel")
    ax.set_yticks(range(8))
    ax.set_ylim(-0.5, 7.5)
    ax.grid(True, alpha=0.3)
    ax.legend(loc="upper right")
    fig.tight_layout()
    fig.savefig(out_path, bbox_inches="tight")
    plt.close(fig)


def plot_bm_case(actual, transformer, bm, out_path, start_hop):
    pos = np.arange(start_hop, start_hop + len(actual))
    fig, ax = plt.subplots(figsize=(13, 4.5), dpi=160)
    ax.step(pos, actual, where="mid", color="#1565C0", linewidth=2, label="Actual next hops")
    ax.scatter(pos, transformer, color="#2E7D32", s=45, marker="o", label="State transformer")
    ax.scatter(pos, bm, color="#C62828", s=55, marker="x", linewidths=2, label="Berlekamp-Massey")
    ax.set_title("Case Study: Berlekamp-Massey Fails While State Transformer Predicts Correctly")
    ax.set_xlabel("Hop index")
    ax.set_ylabel("Frequency channel")
    ax.set_yticks(range(8))
    ax.set_ylim(-0.5, 7.5)
    ax.grid(True, alpha=0.3)
    ax.legend(loc="upper right")
    fig.tight_layout()
    fig.savefig(out_path, bbox_inches="tight")
    plt.close(fig)


def plot_training_history(history, state_dim, out_path, title):
    if not history:
        print(f"Skipping {out_path}: no training history available.")
        return

    epochs = np.arange(1, len(history["train_loss"]) + 1)
    exact_key = "val_seq_acc" if "val_seq_acc" in history else "val_state_acc"
    exact_label = f"Exact {state_dim}-Bit State Acc"

    fig, axes = plt.subplots(1, 3, figsize=(18, 4.5), dpi=160)
    axes[0].plot(epochs, history["train_loss"], color="blue", linewidth=1.8, label="Train Loss")
    if "val_loss" in history:
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

    axes[2].plot(epochs, history[exact_key], color="red", linewidth=1.8, label=exact_label)
    axes[2].set_title(f"Exact {state_dim}-Bit State Accuracy")
    axes[2].set_xlabel("Epoch")
    axes[2].set_ylabel("Accuracy (%)")
    axes[2].set_ylim(0, 102)
    axes[2].grid(True, alpha=0.3)
    axes[2].legend(loc="lower right")

    fig.suptitle(title, y=1.02)
    fig.tight_layout()
    fig.savefig(out_path, bbox_inches="tight")
    plt.close(fig)


def make_gold_plots(args):
    require_file(args.gold_data, "Gold-code dataset")
    require_file(args.gold_checkpoint, "Gold-code state-transformer checkpoint")

    state_dict, checkpoint = load_checkpoint(args.gold_checkpoint)
    model = StateTransformer().to(DEVICE)
    model.load_state_dict(state_dict, strict=True)
    model.eval()

    history = checkpoint.get("history") or load_history(args.gold_history)
    plot_training_history(
        history,
        state_dim=20,
        out_path=os.path.join(args.out_dir, "report_state_transformer_training.png"),
        title="Gold-Code State Transformer Training",
    )

    x_gold, y_gold, snr_gold, _ = load_signal(args.gold_data, args.gold_signal_index)
    _, pred_gold = predict_gold_signal(model, x_gold, batch_size=args.batch_size)
    actual_gold = y_gold[LOOKBACK:]
    gold_acc = float(np.mean(pred_gold == actual_gold) * 100.0)

    plot_spectrogram(
        x_gold,
        y_gold,
        snr_gold,
        os.path.join(args.out_dir, "report_gold_spectrogram.png"),
    )
    plot_prediction(
        actual_gold,
        pred_gold,
        f"Gold Code Signal Prediction (SNR = {snr_gold:.0f} dB, accuracy = {gold_acc:.1f}%)",
        os.path.join(args.out_dir, "report_gold_prediction.png"),
        start_hop=LOOKBACK,
    )

    x_bm, y_bm, _, _ = load_signal(args.gold_data, args.bm_signal_index)
    cnn_preds, pred_bm_model = predict_gold_signal(model, x_bm, batch_size=args.batch_size)
    observed = cnn_preds[: args.bm_observed_hops].copy()
    if args.corrupt_bm_observation >= 0:
        pos = args.corrupt_bm_observation
        if pos < len(observed):
            observed[pos] = (observed[pos] + 7) % 8
    actual_future = y_bm[args.bm_start_hop : args.bm_start_hop + args.bm_future_hops]
    transformer_future = pred_bm_model[
        args.bm_start_hop - LOOKBACK : args.bm_start_hop - LOOKBACK + args.bm_future_hops
    ]
    bm_future = bm_predict_from_hops(observed, num_future=args.bm_future_hops)
    plot_bm_case(
        actual_future,
        transformer_future,
        bm_future,
        os.path.join(args.out_dir, "report_bm_vs_transformer.png"),
        start_hop=args.bm_start_hop,
    )

    return {
        "gold_plot": {
            "signal_index": args.gold_signal_index,
            "snr": snr_gold,
            "accuracy": gold_acc,
        },
        "bm_case": {
            "signal_index": args.bm_signal_index,
            "bm_acc": float(np.mean(bm_future == actual_future) * 100.0),
            "transformer_acc": float(np.mean(transformer_future == actual_future) * 100.0),
        },
    }


def make_mseq_plots(args):
    require_file(args.mseq_data, "M-sequence dataset")
    require_file(args.mseq_checkpoint, "M-sequence state-transformer checkpoint")

    state_dict, checkpoint = load_checkpoint(args.mseq_checkpoint)
    model = MSequenceStateTransformer(state_dim=10).to(DEVICE)
    model.load_state_dict(state_dict, strict=True)
    model.eval()

    history = checkpoint.get("history") or load_history(args.mseq_history)
    plot_training_history(
        history,
        state_dim=10,
        out_path=os.path.join(args.out_dir, "report_mseq_state_transformer_training.png"),
        title="M-Sequence State Transformer Training",
    )

    x_mseq, y_mseq, snr_mseq, seed_mseq = load_signal(args.mseq_data, args.mseq_signal_index)
    pred_mseq = predict_mseq_signal(model, x_mseq, batch_size=args.batch_size)
    actual_mseq = y_mseq[LOOKBACK:]
    mseq_acc = float(np.mean(pred_mseq == actual_mseq) * 100.0)

    plot_prediction(
        actual_mseq,
        pred_mseq,
        f"M-Sequence State-Transformer Prediction (SNR = {snr_mseq:.0f} dB, accuracy = {mseq_acc:.1f}%)",
        os.path.join(args.out_dir, "report_mseq_prediction.png"),
        start_hop=LOOKBACK,
    )

    return {
        "mseq_plot": {
            "signal_index": args.mseq_signal_index,
            "snr": snr_mseq,
            "seed": seed_mseq,
            "accuracy": mseq_acc,
        }
    }


def parse_args():
    parser = argparse.ArgumentParser(description="Generate result plots for FHSS state-recovery experiments.")
    parser.add_argument("--out-dir", default="results")
    parser.add_argument("--batch-size", type=int, default=64)

    parser.add_argument("--gold-data", default=os.path.join("data", "synthetic", "classification_dataset_stft_random_deg10_python.h5"))
    parser.add_argument("--gold-checkpoint", default=os.path.join("models", "state_tracking_transformer_best.pth"))
    parser.add_argument("--gold-history", default=os.path.join("models", "state_transformer_history.json"))
    parser.add_argument("--gold-signal-index", type=int, default=1)

    parser.add_argument("--bm-signal-index", type=int, default=33)
    parser.add_argument("--bm-observed-hops", type=int, default=50)
    parser.add_argument("--bm-start-hop", type=int, default=50)
    parser.add_argument("--bm-future-hops", type=int, default=20)
    parser.add_argument("--corrupt-bm-observation", type=int, default=1)

    parser.add_argument("--mseq-data", default=os.path.join("data", "synthetic", "m_sequence_state_dataset.h5"))
    parser.add_argument("--mseq-checkpoint", default=os.path.join("models", "m_sequence_state_transformer_best.pth"))
    parser.add_argument("--mseq-history", default="")
    parser.add_argument("--mseq-signal-index", type=int, default=42)

    parser.add_argument("--skip-gold", action="store_true")
    parser.add_argument("--skip-mseq", action="store_true")
    return parser.parse_args()


def main():
    args = parse_args()
    ensure_dir(args.out_dir)

    summary = {}
    if not args.skip_gold:
        summary.update(make_gold_plots(args))
    if not args.skip_mseq:
        summary.update(make_mseq_plots(args))

    with open(os.path.join(args.out_dir, "report_figure_summary.json"), "w", encoding="utf-8") as fp:
        json.dump(summary, fp, indent=2)


if __name__ == "__main__":
    main()
