#!/usr/bin/env python3
"""
Summarize CICADA Smoothing_Retention as an isolated feature across smoothing
scales using PRIVATE FIX-training development-set results.

Input:
  smoothing_grid_ic_level_PRIVATE.csv
produced by evaluate_smoothing_grid_fixtraining.py

For each dataset/candidate, reports:
  - manual signal/noise IC counts
  - median Smoothing_Retention for manual signal vs noise
  - pooled descriptive AUC (higher retention predicts manual signal)
  - mean scan-level AUC
  - signal low-smoothing penalty rate:
        P(Low_Smoothing | manual signal)
  - noise escape rate:
        P(not Low_Smoothing | manual noise)
  - signal high-smoothing rate
  - noise high-smoothing rate

This is DEVELOPMENT-set descriptive analysis. ICs are clustered within scans
and are not treated as independent inferential observations.
"""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path
from collections import defaultdict

import numpy as np


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--ic-level", required=True, type=Path)
    p.add_argument("--output", type=Path, default=None)
    return p.parse_args()


def read_rows(path):
    with path.open(newline="", encoding="utf-8-sig") as f:
        rows = list(csv.DictReader(f))
    if not rows:
        raise RuntimeError(f"No rows in {path}")

    required = {
        "sub_id", "ses_id", "task_name", "candidate", "fwhm_mm",
        "manual_label", "smoothing_retention",
        "high_smoothing", "low_smoothing"
    }
    missing = required - set(rows[0])
    if missing:
        raise RuntimeError(f"Missing columns: {sorted(missing)}")
    return rows


def average_ranks(values):
    values = np.asarray(values, dtype=float)
    order = np.argsort(values, kind="stable")
    ranks = np.empty(len(values), dtype=float)

    i = 0
    while i < len(values):
        j = i + 1
        while j < len(values) and values[order[j]] == values[order[i]]:
            j += 1
        avg = (i + 1 + j) / 2.0
        ranks[order[i:j]] = avg
        i = j
    return ranks


def auc_from_scores(truth, scores):
    truth = np.asarray(truth, dtype=bool)
    scores = np.asarray(scores, dtype=float)

    n_pos = int(np.sum(truth))
    n_neg = int(np.sum(~truth))
    if n_pos == 0 or n_neg == 0:
        return math.nan

    ranks = average_ranks(scores)
    sum_pos = np.sum(ranks[truth])
    u = sum_pos - n_pos * (n_pos + 1) / 2.0
    return float(u / (n_pos * n_neg))


def safe_rate(num, den):
    return float(num / den) if den else math.nan


def summarize(rows):
    truth = np.asarray([int(float(r["manual_label"])) for r in rows], dtype=bool)
    retention = np.asarray(
        [float(r["smoothing_retention"]) for r in rows], dtype=float
    )
    low = np.asarray([int(float(r["low_smoothing"])) for r in rows], dtype=bool)
    high = np.asarray([int(float(r["high_smoothing"])) for r in rows], dtype=bool)

    sig = truth
    noise = ~truth

    pooled_auc = auc_from_scores(truth, retention)

    by_scan = defaultdict(list)
    for r in rows:
        key = (r["sub_id"], r["ses_id"], r["task_name"])
        by_scan[key].append(r)

    scan_aucs = []
    for scan_rows in by_scan.values():
        t = np.asarray(
            [int(float(r["manual_label"])) for r in scan_rows], dtype=bool
        )
        s = np.asarray(
            [float(r["smoothing_retention"]) for r in scan_rows], dtype=float
        )
        auc = auc_from_scores(t, s)
        if np.isfinite(auc):
            scan_aucs.append(auc)

    mean_scan_auc = float(np.mean(scan_aucs)) if scan_aucs else math.nan
    sd_scan_auc = (
        float(np.std(scan_aucs, ddof=1)) if len(scan_aucs) > 1 else 0.0
    )

    return {
        "n_scans": len(by_scan),
        "n_ics": len(rows),
        "manual_signal": int(np.sum(sig)),
        "manual_noise": int(np.sum(noise)),
        "signal_retention_median": float(np.median(retention[sig])),
        "noise_retention_median": float(np.median(retention[noise])),
        "median_gap_signal_minus_noise":
            float(np.median(retention[sig]) - np.median(retention[noise])),
        "pooled_auc": pooled_auc,
        "mean_scan_auc": mean_scan_auc,
        "sd_scan_auc": sd_scan_auc,
        "signal_low_smoothing_penalty_rate":
            safe_rate(np.sum(low & sig), np.sum(sig)),
        "noise_escape_low_smoothing_rate":
            safe_rate(np.sum((~low) & noise), np.sum(noise)),
        "signal_high_smoothing_rate":
            safe_rate(np.sum(high & sig), np.sum(sig)),
        "noise_high_smoothing_rate":
            safe_rate(np.sum(high & noise), np.sum(noise)),
    }


def main():
    args = parse_args()
    rows = read_rows(args.ic_level.expanduser().resolve())

    candidates = {}
    for r in rows:
        cname = r["candidate"]
        candidates[cname] = float(r["fwhm_mm"])

    ordered_candidates = sorted(candidates, key=lambda c: candidates[c])

    out_rows = []

    for dataset in ("rest", "foodpics_run-01", "combined_descriptive"):
        for candidate in ordered_candidates:
            subset = [
                r for r in rows
                if r["candidate"] == candidate
                and (
                    dataset == "combined_descriptive"
                    or r["task_name"] == dataset
                )
            ]
            if not subset:
                continue

            stats = summarize(subset)
            out_rows.append({
                "dataset": dataset,
                "candidate": candidate,
                "fwhm_mm": candidates[candidate],
                **stats,
            })

    print("=" * 118)
    print("CICADA SMOOTHING_RETENTION FEATURE CURVE — FIX-TRAINING DEVELOPMENT SET")
    print("=" * 118)

    for dataset in ("rest", "foodpics_run-01", "combined_descriptive"):
        print(f"\n{dataset}")
        print(
            "FWHM    SigMed   NoiseMed  Gap      AUC(pool)  "
            "AUC(scan)  SigLow%  NoiseEscape%  SigHigh%  NoiseHigh%"
        )
        for r in [x for x in out_rows if x["dataset"] == dataset]:
            print(
                f"{r['fwhm_mm']:6.2f}  "
                f"{r['signal_retention_median']:.4f}   "
                f"{r['noise_retention_median']:.4f}    "
                f"{r['median_gap_signal_minus_noise']:+.4f}   "
                f"{r['pooled_auc']:.4f}     "
                f"{r['mean_scan_auc']:.4f}     "
                f"{100*r['signal_low_smoothing_penalty_rate']:6.2f}   "
                f"{100*r['noise_escape_low_smoothing_rate']:10.2f}   "
                f"{100*r['signal_high_smoothing_rate']:7.2f}   "
                f"{100*r['noise_high_smoothing_rate']:9.2f}"
            )

    print(
        "\nInterpretation:"
        "\n  Higher AUC / larger positive median gap = better isolated retention separation."
        "\n  Lower SigLow% = fewer manual signal ICs penalized as low-smoothing."
        "\n  Lower NoiseEscape% = fewer manual noise ICs escape the low-smoothing penalty."
        "\n  A useful operating region should balance both error modes rather than optimize"
        "\n  one metric in isolation."
        "\n\nAUC(pool) is descriptive only because ICs are nested within scans."
        "\nAUC(scan) is the mean of per-scan AUCs and is the more useful robustness view."
    )

    if args.output:
        out = args.output.expanduser().resolve()
        out.parent.mkdir(parents=True, exist_ok=True)
        with out.open("w", newline="", encoding="utf-8") as f:
            fields = list(out_rows[0].keys())
            w = csv.DictWriter(f, fieldnames=fields)
            w.writeheader()
            w.writerows(out_rows)
        print(f"\nPrivate feature-curve summary written to:\n  {out}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
