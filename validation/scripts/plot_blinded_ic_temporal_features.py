#!/usr/bin/env python3
"""
Create blinded time-series and power-spectrum PNGs for CICADA IC review cases.

This script works ONLY from the already-blinded Review_###_timeseries.txt files.
It does not read the unblinding key and therefore does not expose subject,
task, consensus label, historical CICADA label, or 11-mm CICADA label.

For each Review_### case it creates:
  Review_###_timeseries_plot.png
  Review_###_power_spectrum.png

Power spectrum:
  - linear detrending
  - Hann-windowed Welch-style averaged periodogram
  - normalized to unit total power
  - x-axis in Hz using --tr

All outputs remain PRIVATE human-subject-derived review material.
"""

from __future__ import annotations

import argparse
from pathlib import Path
import math

import numpy as np

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument(
        "--review-dir",
        required=True,
        type=Path,
        help="Path containing Review_###_timeseries.txt files.",
    )
    p.add_argument(
        "--tr",
        required=True,
        type=float,
        help="Repetition time in seconds; e.g. 2.0.",
    )
    return p.parse_args()


def load_timeseries(path: Path):
    arr = np.loadtxt(path)

    if arr.ndim == 1:
        return np.asarray(arr, dtype=float)

    if arr.ndim == 2 and arr.shape[1] == 1:
        return np.asarray(arr[:, 0], dtype=float)

    if arr.ndim == 2 and arr.shape[1] == 2:
        first = np.asarray(arr[:, 0], dtype=float)
        second = np.asarray(arr[:, 1], dtype=float)

        diffs = np.diff(first)
        if (
            len(diffs) > 0
            and np.all(np.isfinite(diffs))
            and np.all(diffs > 0)
            and np.std(diffs) < 1e-6 * max(abs(np.mean(diffs)), 1.0)
        ):
            return second

    raise RuntimeError(
        f"Could not safely interpret {path.name}: expected one numeric "
        "time-series column (or a time/index + value pair)."
    )


def linear_detrend(y):
    y = np.asarray(y, dtype=float)
    x = np.arange(len(y), dtype=float)
    design = np.column_stack([np.ones(len(y)), x])
    beta, *_ = np.linalg.lstsq(design, y, rcond=None)
    return y - design @ beta


def zscore(y):
    y = np.asarray(y, dtype=float)
    sd = np.std(y, ddof=1)
    if not np.isfinite(sd) or sd == 0:
        return y - np.mean(y)
    return (y - np.mean(y)) / sd


def choose_segment_length(n):
    # Aim for good readability with typical ~300-volume fMRI runs.
    if n >= 256:
        return 128
    if n >= 128:
        return 64
    if n >= 64:
        return 32
    return max(8, 2 ** int(math.floor(math.log2(max(n, 8)))))


def welch_psd(y, tr):
    y = linear_detrend(y)
    n = len(y)
    nperseg = min(choose_segment_length(n), n)

    if nperseg < 8:
        raise RuntimeError("Time series too short for useful power spectrum.")

    step = max(1, nperseg // 2)
    starts = list(range(0, n - nperseg + 1, step))
    if not starts:
        starts = [0]
        nperseg = n

    window = np.hanning(nperseg)
    window_power = np.sum(window ** 2)

    spectra = []
    for start in starts:
        segment = y[start:start + nperseg]
        segment = linear_detrend(segment)
        fft = np.fft.rfft(segment * window)
        power = (np.abs(fft) ** 2) / window_power
        spectra.append(power)

    psd = np.mean(np.vstack(spectra), axis=0)
    freqs = np.fft.rfftfreq(nperseg, d=tr)

    # Remove DC; it is not useful after detrending and can dominate plotting.
    keep = freqs > 0
    freqs = freqs[keep]
    psd = psd[keep]

    total = np.sum(psd)
    if total > 0:
        psd = psd / total

    return freqs, psd


def save_timeseries_plot(blind_id, y, tr, out_path):
    y_plot = zscore(linear_detrend(y))
    time_seconds = np.arange(len(y_plot)) * tr
    time_minutes = time_seconds / 60.0

    fig, ax = plt.subplots(figsize=(9, 3.5))
    ax.plot(time_minutes, y_plot, linewidth=1.0)
    ax.axhline(0, linewidth=0.7)
    ax.set_title(f"{blind_id} — IC time course")
    ax.set_xlabel("Time (minutes)")
    ax.set_ylabel("Standardized amplitude")
    ax.set_xlim(time_minutes[0], time_minutes[-1] if len(time_minutes) > 1 else 1)
    fig.tight_layout()
    fig.savefig(out_path, dpi=180, bbox_inches="tight")
    plt.close(fig)


def save_power_plot(blind_id, y, tr, out_path):
    freqs, psd = welch_psd(y, tr)

    fig, ax = plt.subplots(figsize=(9, 3.5))
    ax.plot(freqs, psd, linewidth=1.2)
    ax.set_title(f"{blind_id} — normalized temporal power spectrum")
    ax.set_xlabel("Frequency (Hz)")
    ax.set_ylabel("Fraction of total spectral power")
    ax.set_xlim(0, 1.0 / (2.0 * tr))
    ax.set_ylim(bottom=0)
    fig.tight_layout()
    fig.savefig(out_path, dpi=180, bbox_inches="tight")
    plt.close(fig)


def main():
    args = parse_args()
    review_dir = args.review_dir.expanduser().resolve()

    if args.tr <= 0:
        raise RuntimeError("--tr must be > 0 seconds.")

    files = sorted(review_dir.glob("Review_*_timeseries.txt"))
    if not files:
        raise RuntimeError(
            f"No Review_*_timeseries.txt files found in {review_dir}"
        )

    print("=" * 72)
    print("CICADA BLINDED TEMPORAL REVIEW PLOTS")
    print("=" * 72)
    print(f"Review cases: {len(files)}")
    print(f"TR: {args.tr:g} s")
    print(f"Nyquist frequency: {1/(2*args.tr):.4f} Hz")
    print()

    for path in files:
        blind_id = path.name.replace("_timeseries.txt", "")
        y = load_timeseries(path)

        ts_png = review_dir / f"{blind_id}_timeseries_plot.png"
        ps_png = review_dir / f"{blind_id}_power_spectrum.png"

        save_timeseries_plot(blind_id, y, args.tr, ts_png)
        save_power_plot(blind_id, y, args.tr, ps_png)

        print(
            f"{blind_id}: {len(y)} time points -> "
            f"{ts_png.name}, {ps_png.name}"
        )

    print(
        "\nThese plots remain blinded: no consensus or CICADA classification "
        "information was read or displayed."
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
