#!/usr/bin/env python3
"""
CICADA smoothing-only downstream reclassification in Python.

Purpose
-------
Reconstruct CICADA's downstream automatic IC classification from the small,
language-neutral CSV outputs written by CICADA_2_AutoLabeling.m:

  * feature_auto_vals.csv  (feature_relative_table)
  * IC_auto_checker.csv    (frozen automatic labels)

The script first MUST reproduce the frozen baseline SignalLabel values exactly
using the historical Smoothing_Retention values. Only after that passes does it
replace Smoothing_Retention with a corrected vector and rerun the downstream
classifier.

This makes the smoothing experiment independent of a MATLAB license.

IMPORTANT
---------
Outputs are participant-derived research results and are PRIVATE. By default
the script refuses to write inside a Git repository.
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path
import sys

import numpy as np


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--features", required=True, type=Path,
                   help="Frozen feature_auto_vals.csv")
    p.add_argument("--checker", required=True, type=Path,
                   help="Frozen IC_auto_checker.csv")
    p.add_argument("--corrected-retention", required=True, type=Path,
                   help="Corrected Smoothing_Retention vector, one value per IC")
    p.add_argument("--output-dir", required=True, type=Path)
    p.add_argument("--tolerance", type=int, default=5,
                   help="CICADA tolerance used for the frozen run (default 5)")
    p.add_argument("--allow-output-inside-git", action="store_true")
    return p.parse_args()


def find_git_root(path: Path):
    path = path.resolve()
    for candidate in (path, *path.parents):
        if (candidate / ".git").exists():
            return candidate
    return None


def read_numeric_csv(path: Path):
    with path.open(newline="", encoding="utf-8-sig") as f:
        reader = csv.DictReader(f)
        if not reader.fieldnames:
            raise RuntimeError(f"No header found in {path}")
        names = list(reader.fieldnames)
        cols = {name: [] for name in names}
        for row in reader:
            for name in names:
                val = str(row[name]).strip()
                if val == "":
                    cols[name].append(np.nan)
                else:
                    cols[name].append(float(val))
    return names, {k: np.asarray(v, dtype=float) for k, v in cols.items()}


def read_checker(path: Path):
    with path.open(newline="", encoding="utf-8-sig") as f:
        reader = csv.DictReader(f)
        required = {"PotentialICs", "SignalLabel"}
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise RuntimeError(
                f"{path} missing checker columns: {sorted(missing)}"
            )
        out = {}
        for row in reader:
            ic = int(float(row["PotentialICs"]))
            label = int(float(row["SignalLabel"]))
            if label not in (0, 1):
                raise RuntimeError(f"Nonbinary SignalLabel for IC {ic}")
            if ic in out:
                raise RuntimeError(f"Duplicate IC {ic} in {path}")
            out[ic] = label
    return out


def range_normalize(x):
    x = np.asarray(x, dtype=float)
    lo = np.min(x)
    hi = np.max(x)
    if not np.isfinite(lo) or not np.isfinite(hi):
        raise RuntimeError("Non-finite feature value encountered")
    if hi == lo:
        # A constant feature carries no ranking information.
        return np.zeros_like(x)
    return (x - lo) / (hi - lo)


def matlab_style_kmeans_1d(x, max_iter=100):
    """
    Deterministic 1-D Lloyd k-means initialized exactly as CICADA:
        [min(x), median(x), max(x)]

    Returns zero-based cluster assignments and centroids.

    NumPy argmin resolves exact distance ties toward the first centroid, which
    is deterministic. Empty clusters are an error here because the frozen
    baseline reproduction check should expose any behavior mismatch rather
    than silently guessing MATLAB's EmptyAction behavior.
    """
    x = np.asarray(x, dtype=float).reshape(-1)
    if np.any(~np.isfinite(x)):
        raise RuntimeError("k-means feature contains NaN/Inf")

    c = np.array([np.min(x), np.median(x), np.max(x)], dtype=float)

    labels = None
    for _ in range(max_iter):
        distances = np.abs(x[:, None] - c[None, :])
        new_labels = np.argmin(distances, axis=1)

        new_c = c.copy()
        for k in range(3):
            members = x[new_labels == k]
            if members.size == 0:
                raise RuntimeError(
                    "A k-means cluster became empty. "
                    "Python reconstruction cannot safely assume MATLAB's "
                    "empty-cluster handling."
                )
            new_c[k] = np.mean(members)

        if labels is not None and np.array_equal(new_labels, labels):
            c = new_c
            labels = new_labels
            break

        if np.allclose(new_c, c, rtol=0, atol=0):
            c = new_c
            labels = new_labels
            break

        c = new_c
        labels = new_labels
    else:
        raise RuntimeError("1-D k-means failed to converge within 100 iterations")

    # Reassign once to the final centroids, matching the converged partition.
    labels = np.argmin(np.abs(x[:, None] - c[None, :]), axis=1)
    return labels, c


def build_classification(features):
    """
    Reconstruct CICADA classification_table from feature_relative_table.
    """
    if "ICs" not in features:
        raise RuntimeError("feature_auto_vals.csv is missing ICs")

    n = len(features["ICs"])
    classification = {}

    for name, values in features.items():
        if name == "ICs":
            continue
        if len(values) != n:
            raise RuntimeError(f"Feature length mismatch: {name}")

        labels, centers = matlab_style_kmeans_1d(values)
        max_cluster = int(np.argmax(centers))
        min_cluster = int(np.argmin(centers))

        classification[f"High_{name}"] = labels == max_cluster
        classification[f"Low_{name}"] = labels == min_cluster

    # Resting-state CICADA creates this alias so downstream code can use the
    # same rule regardless of task availability.
    if "High_besttask_power_overlap" not in classification:
        if "High_general_power_overlap" not in classification:
            raise RuntimeError("Missing High_general_power_overlap")
        classification["High_besttask_power_overlap"] = (
            classification["High_general_power_overlap"].copy()
        )

    # Historical CICADA combines Subepe and WMCSF evidence.
    for required in ("High_Subepe", "High_WMCSF"):
        if required not in classification:
            raise RuntimeError(f"Missing {required}")
    classification["High_Subepe"] = (
        classification["High_Subepe"] |
        classification["High_WMCSF"]
    )

    # Historical Spikiness gate: only high-kmeans spikiness is bad if raw
    # Spikiness itself exceeds 5.
    if "High_Spikiness" not in classification or "Spikiness" not in features:
        raise RuntimeError("Missing Spikiness feature/classification")
    classification["High_Spikiness"] = (
        classification["High_Spikiness"] &
        (features["Spikiness"] > 5)
    )

    return classification


def classify(features, tolerance):
    classification = build_classification(features)

    best = [
        "High_GM",
        "High_Smoothing_Retention",
        "High_best_power_overlap_norm",
    ]
    worst = [
        "Low_GM",
        "Low_best_power_overlap_norm",
    ]
    bad_region = [
        "Low_GM",
        "High_Edge",
        "High_Subepe",
        "High_CSF",
        "High_Suscept",
        "High_OutbrainOnly",
    ]
    bad = [
        "Low_GM",
        "High_Edge",
        "High_Subepe",
        "High_CSF",
        "High_Suscept",
        "High_OutbrainOnly",
        "High_Outbrain",
        "High_Highfreq",
        "High_Spikiness",
        "Low_best_power_overlap_norm",
        "Low_Smoothing_Retention",
        "High_DVARS_Corr",
        "High_FD_Corr",
    ]

    needed = set(best + worst + bad_region + bad + ["High_besttask_power_overlap"])
    missing = needed - set(classification)
    if missing:
        raise RuntimeError(
            "Cannot reconstruct CICADA decision rules; missing flags: "
            + ", ".join(sorted(missing))
        )

    score = (
        range_normalize(features["Smoothing_Retention"])
        * range_normalize(features["GM"]) ** 2
        * range_normalize(features["best_power_overlap_norm"])
    )

    # Stable descending ordering: for exact ties retain original IC order,
    # matching modern MATLAB stable sort behavior.
    signal_idx0 = np.argsort(-score, kind="stable")
    signal_vals = score[signal_idx0]

    g_end = int(np.sum(signal_vals > np.mean(signal_vals)) + 1)
    stop_num = int(tolerance)

    signal_decision = []
    g0 = 0

    while stop_num > 0 and g0 < (g_end - 1):
        idx = int(signal_idx0[g0])

        high_gm = bool(classification["High_GM"][idx])
        high_best = bool(classification["High_best_power_overlap_norm"][idx])
        high_besttask = bool(classification["High_besttask_power_overlap"][idx])

        worst_sum = sum(bool(classification[k][idx]) for k in worst)
        bad_region_sum = sum(bool(classification[k][idx]) for k in bad_region)
        bad_sum = sum(bool(classification[k][idx]) for k in bad)

        decision = False

        if (high_gm or high_best or high_besttask) and worst_sum == 0:
            if high_gm or bad_region_sum == 0:
                if high_gm and (
                    high_best
                    or bool(classification["High_Smoothing_Retention"][idx])
                ):
                    stop_num += 1
                    decision = True
                elif bad_sum == bad_region_sum:
                    stop_num += 1
                    decision = True
                else:
                    stop_num -= 1
            else:
                stop_num -= 1
        else:
            stop_num -= 1

        signal_decision.append(decision)

        if stop_num > tolerance:
            stop_num = tolerance

        g0 += 1

    potential_idx0 = signal_idx0[:g0]

    n = len(score)
    signal = np.zeros(n, dtype=bool)
    if g0:
        signal[potential_idx0] = np.asarray(signal_decision, dtype=bool)

    # Historical CICADA failsafe exactly as currently coded.
    if np.sum(signal) == 0:
        if n >= 1:
            signal[0] = True
        if n >= 2:
            signal[1] = True

    if np.sum(signal) < 2:
        if n >= 1:
            signal[0] = True
        if np.sum(signal) < 2 and n >= 2:
            signal[1] = True

    highnoise = np.zeros(n, dtype=bool)
    highsig = np.zeros(n, dtype=bool)
    for i in range(n):
        highnoise[i] = any(bool(classification[k][i]) for k in bad)
        highsig[i] = all(bool(classification[k][i]) for k in best)

    return {
        "signal": signal,
        "score": score,
        "signal_order0": signal_idx0,
        "classification": classification,
        "highnoise": highnoise,
        "highsig": highsig,
    }


def write_csv(path, rows, fieldnames):
    with path.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        w.writerows(rows)


def main():
    args = parse_args()

    args.features = args.features.expanduser().resolve()
    args.checker = args.checker.expanduser().resolve()
    args.corrected_retention = args.corrected_retention.expanduser().resolve()
    args.output_dir = args.output_dir.expanduser().resolve()

    git_root = find_git_root(args.output_dir)
    if git_root is not None and not args.allow_output_inside_git:
        raise RuntimeError(
            "\nPRIVACY SAFETY STOP\n"
            "Participant-derived smoothing-validation output was requested "
            "inside a Git repository.\n\n"
            f"Requested output: {args.output_dir}\n"
            f"Git repository:    {git_root}\n"
        )

    args.output_dir.mkdir(parents=True, exist_ok=True)

    names, features = read_numeric_csv(args.features)
    frozen_checker = read_checker(args.checker)

    ics = features.get("ICs")
    if ics is None:
        raise RuntimeError("No ICs column in feature file")
    ics = np.asarray(np.rint(ics), dtype=int)

    if not np.array_equal(ics, np.arange(1, len(ics) + 1)):
        raise RuntimeError(
            "Expected feature_auto_vals.csv ICs to be sequential 1..N"
        )

    corrected = np.loadtxt(args.corrected_retention, dtype=float).reshape(-1)
    if len(corrected) != len(ics):
        raise RuntimeError(
            f"Corrected retention length {len(corrected)} != {len(ics)} ICs"
        )

    if "Smoothing_Retention" not in features:
        raise RuntimeError("Feature file missing Smoothing_Retention")

    # A. Exact baseline reproduction.
    baseline = classify(features, args.tolerance)

    frozen_labels = np.array(
        [frozen_checker.get(int(ic), -1) for ic in ics],
        dtype=int,
    )
    if np.any(frozen_labels < 0):
        missing = ics[frozen_labels < 0]
        raise RuntimeError(
            "Frozen checker missing ICs: " + str(missing.tolist())
        )

    baseline_match = np.array_equal(
        baseline["signal"].astype(int),
        frozen_labels,
    )

    print("=" * 72)
    print("CICADA PYTHON SMOOTHING-ONLY RECLASSIFICATION")
    print("=" * 72)
    print(f"ICs: {len(ics)}")
    print(f"Tolerance: {args.tolerance}")
    print(f"Baseline labels reproduced exactly: {baseline_match}")

    if not baseline_match:
        bad_ics = ics[
            baseline["signal"].astype(int) != frozen_labels
        ]
        print("Mismatching ICs:", bad_ics.tolist())
        print(
            "\nSTOP: Python reconstruction does not exactly reproduce the "
            "frozen CICADA baseline. Do not interpret corrected-smoothing "
            "results."
        )
        return 2

    # B. Replace ONLY smoothing retention.
    corrected_features = {
        k: np.array(v, copy=True) for k, v in features.items()
    }
    corrected_features["Smoothing_Retention"] = corrected

    changed = classify(corrected_features, args.tolerance)

    old_signal = baseline["signal"]
    new_signal = changed["signal"]

    old_high = baseline["classification"]["High_Smoothing_Retention"]
    new_high = changed["classification"]["High_Smoothing_Retention"]
    old_low = baseline["classification"]["Low_Smoothing_Retention"]
    new_low = changed["classification"]["Low_Smoothing_Retention"]

    flips = old_signal != new_signal
    noise_to_signal = (~old_signal) & new_signal
    signal_to_noise = old_signal & (~new_signal)

    print("\nSmoothing cluster membership")
    print("-" * 72)
    print(f"Historical High_Smoothing: {int(np.sum(old_high))}")
    print(f"Corrected  High_Smoothing: {int(np.sum(new_high))}")
    print(
        "High_Smoothing membership changes:",
        int(np.sum(old_high != new_high)),
    )
    print(f"Historical Low_Smoothing:  {int(np.sum(old_low))}")
    print(f"Corrected  Low_Smoothing:  {int(np.sum(new_low))}")
    print(
        "Low_Smoothing membership changes:",
        int(np.sum(old_low != new_low)),
    )

    print("\nFinal automatic classifications")
    print("-" * 72)
    print(f"Baseline signal ICs:  {int(np.sum(old_signal))}")
    print(f"Corrected signal ICs: {int(np.sum(new_signal))}")
    print(f"Total label flips:    {int(np.sum(flips))}")
    print(f"  Noise -> Signal:    {int(np.sum(noise_to_signal))}")
    print(f"  Signal -> Noise:    {int(np.sum(signal_to_noise))}")

    if np.any(flips):
        print("\nFlipped ICs")
        print("IC\tOLD\tNEW\tOLD_RET\tNEW_RET\tDELTA\tOLD_HIGH\tNEW_HIGH\tOLD_LOW\tNEW_LOW")
        for i in np.where(flips)[0]:
            print(
                f"{ics[i]}\t"
                f"{int(old_signal[i])}\t"
                f"{int(new_signal[i])}\t"
                f"{features['Smoothing_Retention'][i]:.6f}\t"
                f"{corrected[i]:.6f}\t"
                f"{corrected[i]-features['Smoothing_Retention'][i]:+.6f}\t"
                f"{int(old_high[i])}\t"
                f"{int(new_high[i])}\t"
                f"{int(old_low[i])}\t"
                f"{int(new_low[i])}"
            )
    else:
        print("\nNo final auto-label flips occurred.")

    rows = []
    for i, ic in enumerate(ics):
        rows.append({
            "PotentialICs": int(ic),
            "BaselineAutoSignalLabel": int(old_signal[i]),
            "CorrectedAutoSignalLabel": int(new_signal[i]),
            "LabelFlipped": int(flips[i]),
            "HistoricalSmoothingRetention":
                f"{features['Smoothing_Retention'][i]:.10f}",
            "CorrectedSmoothingRetention": f"{corrected[i]:.10f}",
            "DeltaSmoothingRetention":
                f"{corrected[i]-features['Smoothing_Retention'][i]:.10f}",
            "HistoricalHighSmoothing": int(old_high[i]),
            "CorrectedHighSmoothing": int(new_high[i]),
            "HistoricalLowSmoothing": int(old_low[i]),
            "CorrectedLowSmoothing": int(new_low[i]),
        })

    write_csv(
        args.output_dir / "smoothing_only_auto_label_comparison_PRIVATE.csv",
        rows,
        list(rows[0].keys()),
    )

    summary = [{
        "NumICs": len(ics),
        "HistoricalHighSmoothing": int(np.sum(old_high)),
        "CorrectedHighSmoothing": int(np.sum(new_high)),
        "HighSmoothingMembershipChanges": int(np.sum(old_high != new_high)),
        "HistoricalLowSmoothing": int(np.sum(old_low)),
        "CorrectedLowSmoothing": int(np.sum(new_low)),
        "LowSmoothingMembershipChanges": int(np.sum(old_low != new_low)),
        "BaselineSignalICs": int(np.sum(old_signal)),
        "CorrectedSignalICs": int(np.sum(new_signal)),
        "TotalLabelFlips": int(np.sum(flips)),
        "NoiseToSignal": int(np.sum(noise_to_signal)),
        "SignalToNoise": int(np.sum(signal_to_noise)),
    }]

    write_csv(
        args.output_dir / "smoothing_only_summary_PRIVATE.csv",
        summary,
        list(summary[0].keys()),
    )

    print(f"\nPrivate outputs written to:\n  {args.output_dir}")
    print(
        "\nIMPORTANT: this technical fixture compares automatic labels only. "
        "Do not use held-out manual truth to tune the smoothing decision."
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
