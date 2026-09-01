#!/usr/bin/env python3
"""
Functional numeric regression for CICADA Smoothing_Retention configuration.

This uses a PRIVATE technical fixture/reference at runtime but contains no
participant-derived data itself.

It verifies that the constants currently present in
basescripts/CICADA_1_MasksandICAs.sh reproduce:

  1. Historical CICADA Smoothing_Retention values and frozen auto labels.
  2. The previously validated 11-mm-FWHM Smoothing_Retention values and
     candidate auto labels recorded BEFORE the source-code patch.

Run from CICADA repository root.
"""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path
import re
import shutil
import subprocess

import numpy as np

import cicada_smoothing_only_reclassification as csr


SOURCE = Path("basescripts/CICADA_1_MasksandICAs.sh")


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--features", required=True, type=Path)
    p.add_argument("--checker", required=True, type=Path)
    p.add_argument("--icthresh", required=True, type=Path)
    p.add_argument("--candidate-reference", required=True, type=Path)
    p.add_argument("--sub-id", required=True)
    p.add_argument("--ses-id", required=True, type=int)
    p.add_argument("--task", required=True)
    p.add_argument("--tolerance", type=int, default=5)
    return p.parse_args()


def require_program(name):
    if shutil.which(name) is None:
        raise RuntimeError(f"Required program not found on PATH: {name}")


def run(cmd):
    return subprocess.run(
        [str(x) for x in cmd],
        check=True,
        capture_output=True,
        text=True,
    ).stdout


def extract_sigmas(text):
    revised = re.search(
        r'revised\)\s*.*?CICADA_SMOOTHING_RETENTION_SIGMA_MM="([^"]+)"',
        text,
        flags=re.S,
    )
    historical = re.search(
        r'historical\)\s*.*?CICADA_SMOOTHING_RETENTION_SIGMA_MM="([^"]+)"',
        text,
        flags=re.S,
    )
    if not revised or not historical:
        raise RuntimeError("Could not parse revised/historical sigma values from source.")
    return float(revised.group(1)), float(historical.group(1))


def parse_fslstats(output):
    means, nvox = [], []
    for line in output.splitlines():
        if not line.strip():
            continue
        parts = line.split()
        if len(parts) < 3:
            raise RuntimeError(f"Unexpected fslstats line: {line!r}")
        means.append(float(parts[0]))
        nvox.append(float(parts[1]))
    return np.asarray(means), np.asarray(nvox)


def retention(icthresh, sigma, tmpdir, label):
    abs_img = tmpdir / f"{label}_abs.nii.gz"
    smooth_img = tmpdir / f"{label}_smooth_abs.nii.gz"

    run(["fslmaths", icthresh, "-abs", abs_img])
    run(["fslmaths", icthresh, "-s", f"{sigma:.12f}", "-abs", smooth_img])

    um, un = parse_fslstats(run(["fslstats", "-t", abs_img, "-M", "-V"]))
    sm, sn = parse_fslstats(run(["fslstats", "-t", smooth_img, "-M", "-V"]))

    u = um * un
    s = sm * sn
    if np.any(u == 0):
        raise RuntimeError("Zero denominator encountered.")
    return s / u


def read_candidate_reference(path, sub_id, ses_id, task):
    with path.open(newline="", encoding="utf-8-sig") as f:
        rows = list(csv.DictReader(f))

    selected = [
        r for r in rows
        if str(r["sub_id"]).strip() == str(sub_id).strip()
        and int(float(r["ses_id"])) == ses_id
        and str(r["task_name"]).strip() == task
    ]
    if not selected:
        raise RuntimeError(
            f"No candidate-reference rows for sub-{sub_id} ses-{ses_id:02d} {task}"
        )

    selected.sort(key=lambda r: int(float(r["PotentialICs"])))
    ics = np.asarray([int(float(r["PotentialICs"])) for r in selected], dtype=int)
    revised_ret = np.asarray(
        [float(r["corrected_smoothing_retention"]) for r in selected],
        dtype=float,
    )
    revised_labels = np.asarray(
        [int(float(r["corrected_auto_label"])) for r in selected],
        dtype=bool,
    )
    baseline_labels = np.asarray(
        [int(float(r["baseline_auto_label"])) for r in selected],
        dtype=bool,
    )
    return ics, revised_ret, revised_labels, baseline_labels


def main():
    args = parse_args()

    if not SOURCE.is_file():
        raise RuntimeError(f"Cannot find {SOURCE}; run from CICADA repo root.")

    require_program("fslmaths")
    require_program("fslstats")

    source_text = SOURCE.read_text(encoding="utf-8")
    revised_sigma, historical_sigma = extract_sigmas(source_text)

    _, features = csr.read_numeric_csv(args.features.expanduser().resolve())
    frozen_checker = csr.read_checker(args.checker.expanduser().resolve())

    feature_ics = np.asarray(np.rint(features["ICs"]), dtype=int)
    frozen_labels = np.asarray(
        [frozen_checker[int(ic)] for ic in feature_ics],
        dtype=bool,
    )

    ref_ics, revised_ref_ret, revised_ref_labels, baseline_ref_labels = (
        read_candidate_reference(
            args.candidate_reference.expanduser().resolve(),
            args.sub_id,
            args.ses_id,
            args.task,
        )
    )

    if not np.array_equal(feature_ics, ref_ics):
        raise RuntimeError("IC identities differ between feature and reference files.")
    if not np.array_equal(frozen_labels, baseline_ref_labels):
        raise RuntimeError("Frozen checker labels differ from pre-patch held-out reference.")

    import tempfile
    with tempfile.TemporaryDirectory(prefix="cicada_numeric_regression_") as td:
        td = Path(td)
        hist_ret = retention(
            args.icthresh.expanduser().resolve(),
            historical_sigma,
            td,
            "historical",
        )
        rev_ret = retention(
            args.icthresh.expanduser().resolve(),
            revised_sigma,
            td,
            "revised",
        )

    historical_feature_ref = np.asarray(
        features["Smoothing_Retention"],
        dtype=float,
    )

    hist_max_abs = float(np.max(np.abs(hist_ret - historical_feature_ref)))
    rev_max_abs = float(np.max(np.abs(rev_ret - revised_ref_ret)))

    hist_features = {k: np.array(v, copy=True) for k, v in features.items()}
    hist_features["Smoothing_Retention"] = hist_ret
    hist_result = csr.classify(hist_features, args.tolerance)["signal"]

    rev_features = {k: np.array(v, copy=True) for k, v in features.items()}
    rev_features["Smoothing_Retention"] = rev_ret
    rev_result = csr.classify(rev_features, args.tolerance)["signal"]

    hist_labels_exact = np.array_equal(hist_result, frozen_labels)
    rev_labels_exact = np.array_equal(rev_result, revised_ref_labels)

    # Candidate reference CSV stores retention rounded to 10 decimal places.
    retention_atol = 2e-9

    print("=" * 82)
    print("CICADA SMOOTHING_RETENTION FUNCTIONAL NUMERIC REGRESSION")
    print("=" * 82)
    print(f"Source revised sigma:    {revised_sigma:.12f} mm")
    print(f"Source historical sigma: {historical_sigma:.12f} mm")
    print(f"ICs: {len(feature_ics)}")
    print()
    print("Historical mode")
    print(f"  max |retention - frozen feature|: {hist_max_abs:.12g}")
    print(f"  frozen auto labels reproduced:   {hist_labels_exact}")
    print()
    print("Revised 11-mm-FWHM mode")
    print(f"  max |retention - prepatch validated reference|: {rev_max_abs:.12g}")
    print(f"  prepatch candidate labels reproduced:           {rev_labels_exact}")
    print()

    ok_hist_ret = hist_max_abs <= retention_atol
    ok_rev_ret = rev_max_abs <= retention_atol

    print(f"Historical numeric regression PASS: {ok_hist_ret and hist_labels_exact}")
    print(f"Revised numeric regression PASS:    {ok_rev_ret and rev_labels_exact}")

    if not (ok_hist_ret and hist_labels_exact and ok_rev_ret and rev_labels_exact):
        raise SystemExit(2)

    print("\nPASS: source configuration reproduces both pre-patch validated behaviors.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
