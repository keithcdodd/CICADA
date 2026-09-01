#!/usr/bin/env python3
"""
CICADA FIX-training smoothing-scale development curve.

Evaluates the historical CICADA smoothness implementation and prespecified
candidate Gaussian smoothing scales on ONLY the historical FIX-training
development set.

Default curve:
    6, 8, 10, 12 mm FWHM
plus the historical CICADA baseline:
    sigma = 6 mm = FWHM 14.12892027 mm

For every scan, the script:
  1. Reconstructs frozen CICADA labels from feature_auto_vals.csv.
  2. Requires exact agreement with IC_auto_checker.csv.
  3. Uses KD individual manual labels as the DEVELOPMENT reference.
  4. Recomputes Smoothing_Retention for each candidate FWHM from
     melodic/ICthresh_zstat.nii.gz.
  5. Replaces ONLY Smoothing_Retention and reruns downstream CICADA logic.
  6. Summarizes performance at the scan level and IC level.

The 30+30 held-out consensus cohort is intentionally not supported here.

All participant-derived outputs are PRIVATE and the script refuses, by default,
to write inside a Git repository.
"""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile

import numpy as np

import cicada_smoothing_only_reclassification as csr


FWHM_PER_SIGMA = 2.354820045
HISTORICAL_SIGMA_MM = 6.0
HISTORICAL_FWHM_MM = HISTORICAL_SIGMA_MM * FWHM_PER_SIGMA


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--master", required=True, type=Path)
    p.add_argument("--data-root", required=True, type=Path)
    p.add_argument("--output-dir", required=True, type=Path)
    p.add_argument("--tolerance", type=int, default=5)
    p.add_argument(
        "--fwhm",
        nargs="+",
        type=float,
        default=[6.0, 8.0, 10.0, 12.0],
        help=(
            "Candidate FWHM values in mm. Historical sigma=6 mm baseline "
            "is always included separately. Default: 6 8 10 12"
        ),
    )
    p.add_argument("--allow-output-inside-git", action="store_true")
    return p.parse_args()


def require_program(name):
    path = shutil.which(name)
    if path is None:
        raise RuntimeError(f"Required program not found on PATH: {name}")
    return path


def run(cmd):
    proc = subprocess.run(
        [str(x) for x in cmd],
        check=True,
        capture_output=True,
        text=True,
    )
    return proc.stdout


def read_master(path):
    with path.open(newline="", encoding="utf-8-sig") as f:
        rows = list(csv.DictReader(f))
    if not rows:
        raise RuntimeError(f"No rows found in master CSV: {path}")
    required = {"sub_id", "ses_id", "task_name", "role", "exclude"}
    missing = required - set(rows[0])
    if missing:
        raise RuntimeError(f"Master CSV missing columns: {sorted(missing)}")
    return rows


def gold_file(scan_dir):
    explicit = scan_dir / "IC_manual_kd_checker.csv"
    legacy = scan_dir / "IC_manual_checker.csv"
    if explicit.is_file():
        return explicit, "KD_individual_explicit_filename"
    if legacy.is_file():
        return legacy, "KD_individual_legacy_default_filename"
    raise FileNotFoundError(
        f"No KD manual label file found in {scan_dir}"
    )


def parse_fslstats_mean_v(output):
    means = []
    nvox = []
    for line in output.splitlines():
        line = line.strip()
        if not line:
            continue
        parts = line.split()
        if len(parts) < 3:
            raise RuntimeError(
                f"Unexpected fslstats -M -V output line: {line!r}"
            )
        means.append(float(parts[0]))
        nvox.append(float(parts[1]))
    if not means:
        raise RuntimeError("No fslstats values parsed")
    return np.asarray(means), np.asarray(nvox)


def sums_from_image(image):
    mean, nvox = parse_fslstats_mean_v(
        run(["fslstats", "-t", image, "-M", "-V"])
    )
    return mean * nvox


def compute_candidate_retentions(icthresh, fwhm_values, tmpdir):
    """
    Compute the common unsmoothed denominator once, then each candidate.
    """
    abs_img = tmpdir / "abs_icthresh.nii.gz"
    run(["fslmaths", icthresh, "-abs", abs_img])
    denominator = sums_from_image(abs_img)

    if np.any(denominator == 0):
        raise RuntimeError("Zero unsmoothed IC sum encountered")

    out = {}

    for fwhm in fwhm_values:
        sigma = fwhm / FWHM_PER_SIGMA
        smoothed = tmpdir / f"smoothed_fwhm_{fwhm:g}.nii.gz"

        run([
            "fslmaths",
            icthresh,
            "-s",
            f"{sigma:.12f}",
            "-abs",
            smoothed,
        ])

        numerator = sums_from_image(smoothed)

        if len(numerator) != len(denominator):
            raise RuntimeError(
                f"Volume-count mismatch at FWHM={fwhm:g} mm"
            )

        out[fwhm] = numerator / denominator

    return out


def binary_metrics(truth, pred):
    truth = np.asarray(truth, dtype=bool)
    pred = np.asarray(pred, dtype=bool)

    tp = int(np.sum(truth & pred))
    fn = int(np.sum(truth & ~pred))
    fp = int(np.sum(~truth & pred))
    tn = int(np.sum(~truth & ~pred))

    def div(a, b):
        return a / b if b else math.nan

    recall = div(tp, tp + fn)
    ppv = div(tp, tp + fp)
    if (
        np.isfinite(recall)
        and np.isfinite(ppv)
        and (recall + ppv) > 0
    ):
        f1 = 2 * recall * ppv / (recall + ppv)
    else:
        f1 = math.nan

    return {
        "n_ics": tp + fn + fp + tn,
        "TP": tp,
        "FN": fn,
        "FP": fp,
        "TN": tn,
        "signal_sensitivity": recall,
        "signal_ppv": ppv,
        "signal_f1": f1,
        "noise_sensitivity": div(tn, tn + fp),
        "npv": div(tn, tn + fn),
        "accuracy": div(tp + tn, tp + tn + fp + fn),
    }


def mean_sd(values):
    vals = np.asarray(
        [v for v in values if np.isfinite(v)],
        dtype=float,
    )
    if vals.size == 0:
        return math.nan, math.nan
    mean = float(np.mean(vals))
    sd = float(np.std(vals, ddof=1)) if vals.size > 1 else 0.0
    return mean, sd


def write_csv(path, rows, fields):
    with path.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)


def candidate_name(fwhm):
    if fwhm is None:
        return "historical_sigma6"
    return f"fwhm_{fwhm:g}mm"


def main():
    args = parse_args()

    args.master = args.master.expanduser().resolve()
    args.data_root = args.data_root.expanduser().resolve()
    args.output_dir = args.output_dir.expanduser().resolve()

    git_root = csr.find_git_root(args.output_dir)
    if git_root is not None and not args.allow_output_inside_git:
        raise RuntimeError(
            "\nPRIVACY SAFETY STOP\n"
            "Participant-derived smoothing-grid output was requested inside "
            "a Git repository.\n\n"
            f"Requested output: {args.output_dir}\n"
            f"Git repository:    {git_root}\n"
        )

    require_program("fslmaths")
    require_program("fslstats")

    fwhm_values = sorted(set(float(x) for x in args.fwhm))
    if any(x <= 0 for x in fwhm_values):
        raise RuntimeError("All candidate FWHM values must be > 0")

    rows = read_master(args.master)
    selected = [
        r for r in rows
        if r.get("role", "").strip() == "fixtraining"
        and r.get("exclude", "0").strip() != "1"
        and r.get("task_name", "").strip()
        in {"rest", "foodpics_run-01"}
    ]

    if len(selected) != 20:
        raise RuntimeError(
            f"Expected 20 FIX-training scans; found {len(selected)}"
        )

    args.output_dir.mkdir(parents=True, exist_ok=True)

    print("=" * 88)
    print("CICADA FIX-TRAINING SMOOTHING-SCALE DEVELOPMENT CURVE")
    print("=" * 88)
    print("Historical baseline:")
    print(
        f"  sigma={HISTORICAL_SIGMA_MM:.3f} mm "
        f"(FWHM={HISTORICAL_FWHM_MM:.5f} mm)"
    )
    print("Candidate FWHM values:", ", ".join(f"{x:g}" for x in fwhm_values))
    print("Selected scans: 20 (10 rest + 10 foodpics)")
    print()

    scan_rows = []
    ic_rows = []
    baseline_failures = []

    metric_names = [
        "signal_sensitivity",
        "signal_ppv",
        "signal_f1",
        "noise_sensitivity",
        "npv",
        "accuracy",
    ]

    for scan_num, m in enumerate(selected, start=1):
        sub = str(m["sub_id"]).strip()
        ses = int(float(m["ses_id"]))
        task = m["task_name"].strip()

        print(f"[{scan_num:02d}/20] sub-{sub} {task}", flush=True)

        scan_dir = (
            args.data_root
            / f"sub-{sub}"
            / f"ses-{ses:02d}"
            / task
            / "ic_auto_selection"
        )
        run_root = scan_dir.parent

        features_file = scan_dir / "feature_auto_vals.csv"
        checker_file = scan_dir / "IC_auto_checker.csv"
        icthresh = run_root / "melodic" / "ICthresh_zstat.nii.gz"
        manual_file, gold_source = gold_file(scan_dir)

        required = [
            features_file,
            checker_file,
            icthresh,
            manual_file,
        ]
        missing = [str(p) for p in required if not p.is_file()]
        if missing:
            raise RuntimeError(
                f"Missing required files for sub-{sub} {task}:\n"
                + "\n".join(missing)
            )

        _, features = csr.read_numeric_csv(features_file)
        frozen_checker = csr.read_checker(checker_file)
        manual_checker = csr.read_checker(manual_file)

        ics = np.asarray(np.rint(features["ICs"]), dtype=int)
        if not np.array_equal(ics, np.arange(1, len(ics) + 1)):
            raise RuntimeError(
                f"Nonsequential feature ICs for sub-{sub} {task}"
            )
        if set(frozen_checker) != set(ics):
            raise RuntimeError(
                f"Auto checker IC-set mismatch for sub-{sub} {task}"
            )
        if set(manual_checker) != set(ics):
            raise RuntimeError(
                f"Manual checker IC-set mismatch for sub-{sub} {task}"
            )

        frozen_labels = np.asarray(
            [frozen_checker[int(ic)] for ic in ics],
            dtype=bool,
        )
        manual_labels = np.asarray(
            [manual_checker[int(ic)] for ic in ics],
            dtype=bool,
        )

        baseline = csr.classify(features, args.tolerance)
        if not np.array_equal(baseline["signal"], frozen_labels):
            mismatches = ics[
                baseline["signal"] != frozen_labels
            ].tolist()
            baseline_failures.append({
                "sub_id": sub,
                "ses_id": ses,
                "task_name": task,
                "mismatching_ics": ";".join(map(str, mismatches)),
            })
            print(f"  BASELINE FAILURE: {mismatches}")
            continue

        # Candidate retentions from the frozen IC maps.
        with tempfile.TemporaryDirectory(
            prefix="cicada_smoothing_grid_"
        ) as td:
            retentions = compute_candidate_retentions(
                icthresh,
                fwhm_values,
                Path(td),
            )

        candidate_results = {
            None: baseline,
        }

        for fwhm in fwhm_values:
            r = retentions[fwhm]
            if len(r) != len(ics):
                raise RuntimeError(
                    f"Retention length mismatch at {fwhm:g} mm"
                )

            candidate_features = {
                k: np.array(v, copy=True)
                for k, v in features.items()
            }
            candidate_features["Smoothing_Retention"] = r
            candidate_results[fwhm] = csr.classify(
                candidate_features,
                args.tolerance,
            )

        base_signal = baseline["signal"]

        # One scan-level row per candidate, including historical baseline.
        for fwhm in [None] + fwhm_values:
            result = candidate_results[fwhm]
            pred = result["signal"]
            metrics = binary_metrics(manual_labels, pred)

            flips_vs_hist = pred != base_signal
            corrected_errors = (
                (base_signal != manual_labels)
                & (pred == manual_labels)
            )
            introduced_errors = (
                (base_signal == manual_labels)
                & (pred != manual_labels)
            )

            scan_row = {
                "sub_id": sub,
                "ses_id": ses,
                "task_name": task,
                "gold_label_source": gold_source,
                "candidate": candidate_name(fwhm),
                "fwhm_mm":
                    HISTORICAL_FWHM_MM if fwhm is None else fwhm,
                "sigma_mm":
                    HISTORICAL_SIGMA_MM
                    if fwhm is None
                    else fwhm / FWHM_PER_SIGMA,
                "is_historical_baseline": int(fwhm is None),
                "n_ics": len(ics),
                "label_flips_vs_historical":
                    int(np.sum(flips_vs_hist)),
                "manual_errors_corrected_vs_historical":
                    int(np.sum(corrected_errors)),
                "manual_errors_introduced_vs_historical":
                    int(np.sum(introduced_errors)),
            }
            for metric in metric_names:
                scan_row[metric] = metrics[metric]

            scan_rows.append(scan_row)

            if fwhm is None:
                retention_vec = features["Smoothing_Retention"]
            else:
                retention_vec = retentions[fwhm]

            high_s = result["classification"]["High_Smoothing_Retention"]
            low_s = result["classification"]["Low_Smoothing_Retention"]

            for i, ic in enumerate(ics):
                ic_rows.append({
                    "sub_id": sub,
                    "ses_id": ses,
                    "task_name": task,
                    "gold_label_source": gold_source,
                    "candidate": candidate_name(fwhm),
                    "fwhm_mm":
                        HISTORICAL_FWHM_MM if fwhm is None else fwhm,
                    "sigma_mm":
                        HISTORICAL_SIGMA_MM
                        if fwhm is None
                        else fwhm / FWHM_PER_SIGMA,
                    "PotentialICs": int(ic),
                    "manual_label": int(manual_labels[i]),
                    "historical_auto_label": int(base_signal[i]),
                    "candidate_auto_label": int(pred[i]),
                    "label_flipped_vs_historical":
                        int(pred[i] != base_signal[i]),
                    "manual_error_corrected_vs_historical":
                        int(
                            (base_signal[i] != manual_labels[i])
                            and (pred[i] == manual_labels[i])
                        ),
                    "manual_error_introduced_vs_historical":
                        int(
                            (base_signal[i] == manual_labels[i])
                            and (pred[i] != manual_labels[i])
                        ),
                    "smoothing_retention":
                        f"{retention_vec[i]:.10f}",
                    "high_smoothing": int(high_s[i]),
                    "low_smoothing": int(low_s[i]),
                })

        candidate_flip_counts = []
        for fwhm in fwhm_values:
            flips = int(np.sum(
                candidate_results[fwhm]["signal"] != base_signal
            ))
            candidate_flip_counts.append(f"{fwhm:g}mm:{flips}")
        print(
            "  baseline exact; flips vs historical = "
            + ", ".join(candidate_flip_counts)
        )

    if baseline_failures:
        failure_path = (
            args.output_dir
            / "BASELINE_REPRODUCTION_FAILURES_PRIVATE.csv"
        )
        write_csv(
            failure_path,
            baseline_failures,
            ["sub_id", "ses_id", "task_name", "mismatching_ics"],
        )
        print("\nSTOP: at least one historical baseline failed exact reproduction.")
        print(f"See: {failure_path}")
        return 2

    # Cohort guard.
    baseline_scans = [
        r for r in scan_rows
        if r["is_historical_baseline"] == 1
    ]
    if (
        len([r for r in baseline_scans if r["task_name"] == "rest"]) != 10
        or len([
            r for r in baseline_scans
            if r["task_name"] == "foodpics_run-01"
        ]) != 10
    ):
        raise RuntimeError("Did not recover exact 10 + 10 development scans")

    # Mean-scan summary (primary descriptive view).
    summary_rows = []

    candidate_order = [None] + fwhm_values
    for task_group in ["rest", "foodpics_run-01", "combined_descriptive"]:
        for fwhm in candidate_order:
            cname = candidate_name(fwhm)
            subset = [
                r for r in scan_rows
                if r["candidate"] == cname
                and (
                    task_group == "combined_descriptive"
                    or r["task_name"] == task_group
                )
            ]

            row = {
                "dataset": task_group,
                "candidate": cname,
                "fwhm_mm":
                    HISTORICAL_FWHM_MM if fwhm is None else fwhm,
                "sigma_mm":
                    HISTORICAL_SIGMA_MM
                    if fwhm is None
                    else fwhm / FWHM_PER_SIGMA,
                "is_historical_baseline": int(fwhm is None),
                "n_scans": len(subset),
                "total_label_flips_vs_historical":
                    int(sum(r["label_flips_vs_historical"] for r in subset)),
                "manual_errors_corrected_vs_historical":
                    int(sum(
                        r["manual_errors_corrected_vs_historical"]
                        for r in subset
                    )),
                "manual_errors_introduced_vs_historical":
                    int(sum(
                        r["manual_errors_introduced_vs_historical"]
                        for r in subset
                    )),
            }

            for metric in metric_names:
                mean, sd = mean_sd([r[metric] for r in subset])
                row[f"{metric}_mean"] = mean
                row[f"{metric}_sd"] = sd

            summary_rows.append(row)

    # Add deltas to the historical baseline within each task group.
    for task_group in ["rest", "foodpics_run-01", "combined_descriptive"]:
        base = next(
            r for r in summary_rows
            if r["dataset"] == task_group
            and r["is_historical_baseline"] == 1
        )
        for row in summary_rows:
            if row["dataset"] != task_group:
                continue
            for metric in metric_names:
                row[f"delta_{metric}_mean_vs_historical"] = (
                    row[f"{metric}_mean"] - base[f"{metric}_mean"]
                )

    # Private writes.
    write_csv(
        args.output_dir / "smoothing_grid_scan_level_PRIVATE.csv",
        scan_rows,
        list(scan_rows[0].keys()),
    )
    write_csv(
        args.output_dir / "smoothing_grid_ic_level_PRIVATE.csv",
        ic_rows,
        list(ic_rows[0].keys()),
    )
    write_csv(
        args.output_dir / "smoothing_grid_mean_scan_summary_PRIVATE.csv",
        summary_rows,
        list(summary_rows[0].keys()),
    )

    # Aggregate terminal output only; no participant IDs.
    print("\n" + "=" * 110)
    print("MEAN SCAN-LEVEL DEVELOPMENT CURVE")
    print("=" * 110)

    for task_group in ["rest", "foodpics_run-01", "combined_descriptive"]:
        print(f"\n{task_group}")
        print(
            "FWHM(mm)   F1          Sens        NoiseSens   Accuracy    "
            "Flips  Corrected  Introduced"
        )
        for row in [
            x for x in summary_rows if x["dataset"] == task_group
        ]:
            label = (
                f"{row['fwhm_mm']:.2f}*"
                if row["is_historical_baseline"]
                else f"{row['fwhm_mm']:.2f}"
            )
            print(
                f"{label:10s} "
                f"{row['signal_f1_mean']:.4f} "
                f"({row['delta_signal_f1_mean_vs_historical']:+.4f})  "
                f"{row['signal_sensitivity_mean']:.4f} "
                f"({row['delta_signal_sensitivity_mean_vs_historical']:+.4f})  "
                f"{row['noise_sensitivity_mean']:.4f} "
                f"({row['delta_noise_sensitivity_mean_vs_historical']:+.4f})  "
                f"{row['accuracy_mean']:.4f} "
                f"({row['delta_accuracy_mean_vs_historical']:+.4f})  "
                f"{row['total_label_flips_vs_historical']:5d}  "
                f"{row['manual_errors_corrected_vs_historical']:9d}  "
                f"{row['manual_errors_introduced_vs_historical']:10d}"
            )

    print("\n* Historical CICADA baseline: sigma=6 mm, "
          f"FWHM={HISTORICAL_FWHM_MM:.5f} mm.")
    print("\nAll 20 frozen automatic baselines reproduced exactly: True")
    print(f"Private outputs:\n  {args.output_dir}")
    print(
        "\nDo not select a kernel solely by the single highest metric. "
        "Inspect the shape/stability of the curve, rest/task consistency, "
        "and the direction of changed ICs."
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
