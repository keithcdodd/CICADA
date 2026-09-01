#!/usr/bin/env python3
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
CANDIDATE_FWHM_MM = 11.0
CANDIDATE_SIGMA_MM = CANDIDATE_FWHM_MM / FWHM_PER_SIGMA
HISTORICAL_SIGMA_MM = 6.0
HISTORICAL_FWHM_MM = HISTORICAL_SIGMA_MM * FWHM_PER_SIGMA


def parse_args():
    p = argparse.ArgumentParser(
        description="Locked held-out CICADA validation: historical smoothing vs frozen 11-mm FWHM candidate."
    )
    p.add_argument("--master", required=True, type=Path)
    p.add_argument("--data-root", required=True, type=Path)
    p.add_argument("--output-dir", required=True, type=Path)
    p.add_argument("--tolerance", type=int, default=5)
    p.add_argument("--allow-output-inside-git", action="store_true")
    return p.parse_args()


def run(cmd):
    return subprocess.run(
        [str(x) for x in cmd], check=True, capture_output=True, text=True
    ).stdout


def require_program(name):
    if shutil.which(name) is None:
        raise RuntimeError(f"Required program not found on PATH: {name}")


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


def consensus_gold_file(scan_dir):
    p = scan_dir / "IC_manual_checker.csv"
    if not p.is_file():
        raise FileNotFoundError(f"Held-out consensus label file not found: {p}")
    return p


def parse_fslstats_mean_v(output):
    means, nvox = [], []
    for line in output.splitlines():
        line = line.strip()
        if not line:
            continue
        parts = line.split()
        if len(parts) < 3:
            raise RuntimeError(f"Unexpected fslstats output: {line!r}")
        means.append(float(parts[0]))
        nvox.append(float(parts[1]))
    if not means:
        raise RuntimeError("No fslstats values parsed")
    return np.asarray(means), np.asarray(nvox)


def candidate_retention_from_icthresh(icthresh, tmpdir):
    abs_img = tmpdir / "abs_icthresh.nii.gz"
    smoothed = tmpdir / "smoothed_11mmFWHM_abs.nii.gz"
    run(["fslmaths", icthresh, "-abs", abs_img])
    run([
        "fslmaths", icthresh, "-s", f"{CANDIDATE_SIGMA_MM:.12f}",
        "-abs", smoothed
    ])
    u_mean, u_nvox = parse_fslstats_mean_v(
        run(["fslstats", "-t", abs_img, "-M", "-V"])
    )
    s_mean, s_nvox = parse_fslstats_mean_v(
        run(["fslstats", "-t", smoothed, "-M", "-V"])
    )
    u_sum = u_mean * u_nvox
    s_sum = s_mean * s_nvox
    if np.any(u_sum == 0):
        raise RuntimeError("Zero unsmoothed IC sum encountered")
    if len(u_sum) != len(s_sum):
        raise RuntimeError("Volume-count mismatch")
    return s_sum / u_sum


def binary_metrics(truth, pred):
    truth = np.asarray(truth, dtype=bool)
    pred = np.asarray(pred, dtype=bool)
    tp = int(np.sum(truth & pred))
    fn = int(np.sum(truth & ~pred))
    fp = int(np.sum(~truth & pred))
    tn = int(np.sum(~truth & ~pred))

    def div(a, b):
        return a / b if b else math.nan

    sens = div(tp, tp + fn)
    ppv = div(tp, tp + fp)
    f1 = (
        2 * sens * ppv / (sens + ppv)
        if np.isfinite(sens) and np.isfinite(ppv) and (sens + ppv) > 0
        else math.nan
    )
    return {
        "signal_sensitivity": sens,
        "signal_ppv": ppv,
        "signal_f1": f1,
        "noise_sensitivity": div(tn, tn + fp),
        "npv": div(tn, tn + fn),
        "accuracy": div(tp + tn, tp + tn + fp + fn),
    }


def mean_sd(values):
    vals = np.asarray([v for v in values if np.isfinite(v)], dtype=float)
    if vals.size == 0:
        return math.nan, math.nan
    return float(np.mean(vals)), (
        float(np.std(vals, ddof=1)) if vals.size > 1 else 0.0
    )


def write_csv(path, rows):
    with path.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)


def main():
    args = parse_args()
    args.master = args.master.expanduser().resolve()
    args.data_root = args.data_root.expanduser().resolve()
    args.output_dir = args.output_dir.expanduser().resolve()

    git_root = csr.find_git_root(args.output_dir)
    if git_root is not None and not args.allow_output_inside_git:
        raise RuntimeError(
            "\nPRIVACY SAFETY STOP\n"
            "Held-out participant-derived output was requested inside Git.\n"
            f"Requested output: {args.output_dir}\n"
            f"Git repository:    {git_root}\n"
        )

    require_program("fslmaths")
    require_program("fslstats")

    rows = read_master(args.master)
    selected = [
        r for r in rows
        if r.get("role", "").strip() == "manual"
        and r.get("exclude", "0").strip() != "1"
        and r.get("task_name", "").strip() in {"rest", "foodpics_run-01"}
    ]
    rest = [r for r in selected if r["task_name"].strip() == "rest"]
    task = [r for r in selected if r["task_name"].strip() == "foodpics_run-01"]

    print("=" * 88)
    print("CICADA LOCKED HELD-OUT SMOOTHING VALIDATION")
    print("=" * 88)
    print(f"Historical: sigma=6.000 mm (FWHM={HISTORICAL_FWHM_MM:.5f} mm)")
    print(
        f"Frozen candidate: FWHM=11.00 mm "
        f"(sigma={CANDIDATE_SIGMA_MM:.9f} mm)"
    )
    print(f"Held-out scans selected: {len(selected)} ({len(rest)} rest + {len(task)} foodpics)")
    print("Reference: KD/LS/MM consensus via IC_manual_checker.csv")
    print()

    if len(selected) != 60 or len(rest) != 30 or len(task) != 30:
        raise RuntimeError(
            "LOCKED COHORT STOP: expected exactly 60 held-out scans "
            "(30 rest + 30 foodpics)."
        )

    args.output_dir.mkdir(parents=True, exist_ok=True)

    metric_names = [
        "signal_sensitivity", "signal_ppv", "signal_f1",
        "noise_sensitivity", "npv", "accuracy"
    ]
    scan_rows, ic_rows, failures = [], [], []

    for idx, m in enumerate(selected, start=1):
        sub = str(m["sub_id"]).strip()
        ses = int(float(m["ses_id"]))
        task_name = m["task_name"].strip()
        print(f"[{idx:02d}/60] sub-{sub} {task_name} ...", flush=True)

        scan_dir = (
            args.data_root / f"sub-{sub}" / f"ses-{ses:02d}" /
            task_name / "ic_auto_selection"
        )
        run_root = scan_dir.parent
        features_file = scan_dir / "feature_auto_vals.csv"
        checker_file = scan_dir / "IC_auto_checker.csv"
        manual_file = consensus_gold_file(scan_dir)
        icthresh = run_root / "melodic" / "ICthresh_zstat.nii.gz"

        for p in (features_file, checker_file, manual_file, icthresh):
            if not p.is_file():
                raise FileNotFoundError(f"Missing required file: {p}")

        _, features = csr.read_numeric_csv(features_file)
        frozen_checker = csr.read_checker(checker_file)
        manual_checker = csr.read_checker(manual_file)

        ics = np.asarray(np.rint(features["ICs"]), dtype=int)
        if not np.array_equal(ics, np.arange(1, len(ics) + 1)):
            raise RuntimeError(f"Nonsequential ICs for sub-{sub} {task_name}")
        if set(frozen_checker) != set(ics) or set(manual_checker) != set(ics):
            raise RuntimeError(f"IC-set mismatch for sub-{sub} {task_name}")

        frozen_labels = np.asarray([frozen_checker[int(ic)] for ic in ics], dtype=bool)
        manual_labels = np.asarray([manual_checker[int(ic)] for ic in ics], dtype=bool)

        historical = csr.classify(features, args.tolerance)
        historical_signal = historical["signal"]

        if not np.array_equal(historical_signal, frozen_labels):
            mismatches = ics[historical_signal != frozen_labels].tolist()
            failures.append({
                "sub_id": sub, "ses_id": ses, "task_name": task_name,
                "mismatching_ics": ";".join(map(str, mismatches))
            })
            print(f"  BASELINE REPRODUCTION FAILED: {mismatches}")
            continue

        with tempfile.TemporaryDirectory(prefix="cicada_heldout_11mm_") as td:
            candidate_retention = candidate_retention_from_icthresh(
                icthresh, Path(td)
            )

        candidate_features = {k: np.array(v, copy=True) for k, v in features.items()}
        candidate_features["Smoothing_Retention"] = candidate_retention
        candidate = csr.classify(candidate_features, args.tolerance)
        candidate_signal = candidate["signal"]

        flips = historical_signal != candidate_signal
        hist_err = historical_signal != manual_labels
        cand_err = candidate_signal != manual_labels
        corrected = hist_err & ~cand_err
        introduced = ~hist_err & cand_err

        hm = binary_metrics(manual_labels, historical_signal)
        cm = binary_metrics(manual_labels, candidate_signal)

        scan_row = {
            "sub_id": sub,
            "ses_id": ses,
            "task_name": task_name,
            "gold_label_source": "consensus_KD_LS_MM",
            "n_ics": len(ics),
            "baseline_exact_reproduction": 1,
            "historical_signal_ics": int(np.sum(historical_signal)),
            "candidate_signal_ics": int(np.sum(candidate_signal)),
            "total_label_flips": int(np.sum(flips)),
            "noise_to_signal": int(np.sum(~historical_signal & candidate_signal)),
            "signal_to_noise": int(np.sum(historical_signal & ~candidate_signal)),
            "manual_errors_corrected": int(np.sum(corrected)),
            "manual_errors_introduced": int(np.sum(introduced)),
        }
        for metric in metric_names:
            scan_row[f"historical_{metric}"] = hm[metric]
            scan_row[f"candidate_{metric}"] = cm[metric]
        scan_rows.append(scan_row)

        hist_high = historical["classification"]["High_Smoothing_Retention"]
        cand_high = candidate["classification"]["High_Smoothing_Retention"]
        hist_low = historical["classification"]["Low_Smoothing_Retention"]
        cand_low = candidate["classification"]["Low_Smoothing_Retention"]

        for i, ic in enumerate(ics):
            ic_rows.append({
                "sub_id": sub,
                "ses_id": ses,
                "task_name": task_name,
                "gold_label_source": "consensus_KD_LS_MM",
                "PotentialICs": int(ic),
                "manual_label": int(manual_labels[i]),
                "baseline_auto_label": int(historical_signal[i]),
                "corrected_auto_label": int(candidate_signal[i]),
                "label_flipped": int(flips[i]),
                "manual_error_corrected": int(corrected[i]),
                "manual_error_introduced": int(introduced[i]),
                "historical_smoothing_retention": f"{features['Smoothing_Retention'][i]:.10f}",
                "corrected_smoothing_retention": f"{candidate_retention[i]:.10f}",
                "delta_smoothing_retention": f"{candidate_retention[i]-features['Smoothing_Retention'][i]:.10f}",
                "historical_high_smoothing": int(hist_high[i]),
                "corrected_high_smoothing": int(cand_high[i]),
                "historical_low_smoothing": int(hist_low[i]),
                "corrected_low_smoothing": int(cand_low[i]),
            })

        print(
            f"  baseline exact; flips={int(np.sum(flips))}; "
            f"consensus errors corrected={int(np.sum(corrected))}; "
            f"errors introduced={int(np.sum(introduced))}"
        )

    if failures:
        failure_path = args.output_dir / "BASELINE_REPRODUCTION_FAILURES_PRIVATE.csv"
        write_csv(failure_path, failures)
        print("\nSTOP: at least one held-out baseline failed exact reproduction.")
        print(f"See: {failure_path}")
        return 2

    if len(scan_rows) != 60:
        raise RuntimeError(f"Expected 60 successful scans; got {len(scan_rows)}")

    pooled_rows, macro_rows = [], []

    for dataset, task_filter in [
        ("rest", "rest"),
        ("foodpics_run-01", "foodpics_run-01"),
        ("combined_descriptive", None),
    ]:
        ic_subset = [
            r for r in ic_rows if task_filter is None or r["task_name"] == task_filter
        ]
        scan_subset = [
            r for r in scan_rows if task_filter is None or r["task_name"] == task_filter
        ]

        truth = np.asarray([r["manual_label"] for r in ic_subset], dtype=bool)
        hist_pred = np.asarray([r["baseline_auto_label"] for r in ic_subset], dtype=bool)
        cand_pred = np.asarray([r["corrected_auto_label"] for r in ic_subset], dtype=bool)

        hm = binary_metrics(truth, hist_pred)
        cm = binary_metrics(truth, cand_pred)

        pooled = {
            "dataset": dataset,
            "n_scans": len(scan_subset),
            "n_ics": len(ic_subset),
            "label_flips": int(np.sum(hist_pred != cand_pred)),
            "consensus_errors_corrected": int(np.sum((hist_pred != truth) & (cand_pred == truth))),
            "consensus_errors_introduced": int(np.sum((hist_pred == truth) & (cand_pred != truth))),
        }
        for metric in metric_names:
            pooled[f"historical_{metric}"] = hm[metric]
            pooled[f"candidate_{metric}"] = cm[metric]
            pooled[f"delta_{metric}"] = cm[metric] - hm[metric]
        pooled_rows.append(pooled)

        macro = {
            "dataset": dataset,
            "n_scans": len(scan_subset),
            "total_label_flips": int(sum(r["total_label_flips"] for r in scan_subset)),
            "consensus_errors_corrected": int(sum(r["manual_errors_corrected"] for r in scan_subset)),
            "consensus_errors_introduced": int(sum(r["manual_errors_introduced"] for r in scan_subset)),
        }
        for metric in metric_names:
            hmean, hsd = mean_sd([r[f"historical_{metric}"] for r in scan_subset])
            cmean, csd = mean_sd([r[f"candidate_{metric}"] for r in scan_subset])
            macro[f"historical_{metric}_mean"] = hmean
            macro[f"historical_{metric}_sd"] = hsd
            macro[f"candidate_{metric}_mean"] = cmean
            macro[f"candidate_{metric}_sd"] = csd
            macro[f"delta_{metric}_mean"] = cmean - hmean
        macro_rows.append(macro)

    write_csv(args.output_dir / "heldout_11mm_ic_level_PRIVATE.csv", ic_rows)
    write_csv(args.output_dir / "heldout_11mm_scan_level_PRIVATE.csv", scan_rows)
    write_csv(args.output_dir / "heldout_11mm_pooled_summary_PRIVATE.csv", pooled_rows)
    write_csv(args.output_dir / "heldout_11mm_mean_scan_summary_PRIVATE.csv", macro_rows)

    print("\n" + "=" * 96)
    print("LOCKED HELD-OUT CONSENSUS RESULT")
    print("=" * 96)
    print("All 60 frozen automatic baselines reproduced exactly: True")

    for row in macro_rows:
        print(
            f"\n{row['dataset']}: scans={row['n_scans']}, "
            f"label flips={row['total_label_flips']}, "
            f"consensus errors corrected={row['consensus_errors_corrected']}, "
            f"errors introduced={row['consensus_errors_introduced']}"
        )
        print(
            f"  Accuracy: {row['historical_accuracy_mean']:.4f} -> "
            f"{row['candidate_accuracy_mean']:.4f} "
            f"(delta {row['delta_accuracy_mean']:+.4f})"
        )
        print(
            f"  Signal F1: {row['historical_signal_f1_mean']:.4f} -> "
            f"{row['candidate_signal_f1_mean']:.4f} "
            f"(delta {row['delta_signal_f1_mean']:+.4f})"
        )
        print(
            f"  Signal sensitivity: "
            f"{row['historical_signal_sensitivity_mean']:.4f} -> "
            f"{row['candidate_signal_sensitivity_mean']:.4f} "
            f"(delta {row['delta_signal_sensitivity_mean']:+.4f})"
        )
        print(
            f"  Noise sensitivity: "
            f"{row['historical_noise_sensitivity_mean']:.4f} -> "
            f"{row['candidate_noise_sensitivity_mean']:.4f} "
            f"(delta {row['delta_noise_sensitivity_mean']:+.4f})"
        )

    print(f"\nPrivate outputs:\n  {args.output_dir}")
    print(
        "\nThis is the locked held-out test. Do not use these results to "
        "select another FWHM. If 11 mm is rejected, retain historical behavior."
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
