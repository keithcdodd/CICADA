#!/usr/bin/env python3
"""
Create a PRIVATE, BLINDED visual-review packet for final-label IC flips from
the locked held-out CICADA smoothing validation.

Input:
  heldout_11mm_ic_level_PRIVATE.csv

The blinded folder contains neutral IDs only:
  Review_001_ICmap.nii.gz
  Review_001_timeseries.txt
  ...
  BLINDED_SCORING.csv
  README_BLINDED.txt

The separate unblinding folder contains:
  PRIVATE_UNBLINDING_KEY.csv

All outputs are PRIVATE human-subject-derived validation material.
Do not commit or upload them to public GitHub.
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path
import random
import secrets
import shutil
import subprocess


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--comparison", required=True, type=Path)
    p.add_argument("--data-root", required=True, type=Path)
    p.add_argument("--output-dir", required=True, type=Path)
    p.add_argument("--seed", type=int, default=None)
    p.add_argument("--allow-output-inside-git", action="store_true")
    return p.parse_args()


def find_git_root(path: Path):
    path = path.resolve()
    for candidate in (path, *path.parents):
        if (candidate / ".git").exists():
            return candidate
    return None


def require_program(name):
    if shutil.which(name) is None:
        raise RuntimeError(f"Required program not found on PATH: {name}")


def write_csv(path, rows, fields):
    with path.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)


def main():
    args = parse_args()
    args.comparison = args.comparison.expanduser().resolve()
    args.data_root = args.data_root.expanduser().resolve()
    args.output_dir = args.output_dir.expanduser().resolve()

    git_root = find_git_root(args.output_dir)
    if git_root is not None and not args.allow_output_inside_git:
        raise RuntimeError(
            "\nPRIVACY SAFETY STOP\n"
            "Participant-derived blinded review output was requested inside Git.\n\n"
            f"Requested output: {args.output_dir}\n"
            f"Git repository:    {git_root}\n"
        )

    require_program("fslroi")

    with args.comparison.open(newline="", encoding="utf-8-sig") as f:
        reader = csv.DictReader(f)
        rows = list(reader)
        fields = set(reader.fieldnames or [])

    required = {
        "sub_id", "ses_id", "task_name", "PotentialICs",
        "manual_label", "baseline_auto_label", "corrected_auto_label",
        "label_flipped", "manual_error_corrected",
        "manual_error_introduced", "historical_smoothing_retention",
        "corrected_smoothing_retention", "delta_smoothing_retention",
        "gold_label_source",
    }
    missing = required - fields
    if missing:
        raise RuntimeError(
            "Comparison CSV missing required columns: "
            + ", ".join(sorted(missing))
        )

    flipped = [r for r in rows if int(float(r["label_flipped"])) == 1]
    if not flipped:
        raise RuntimeError("No final-label-flipped ICs found.")

    seed = args.seed if args.seed is not None else secrets.randbits(63)
    rng = random.Random(seed)
    rng.shuffle(flipped)

    blinded_dir = args.output_dir / "blinded_review"
    key_dir = args.output_dir / "unblinding"

    if blinded_dir.exists() or key_dir.exists():
        raise RuntimeError(
            "Review directories already exist. Choose a new --output-dir."
        )

    blinded_dir.mkdir(parents=True)
    key_dir.mkdir(parents=True)

    scoring_rows = []
    key_rows = []

    for number, row in enumerate(flipped, start=1):
        blind_id = f"Review_{number:03d}"

        sub = str(row["sub_id"]).strip()
        ses = int(float(row["ses_id"]))
        task = str(row["task_name"]).strip()
        ic = int(float(row["PotentialICs"]))

        run_root = (
            args.data_root
            / f"sub-{sub}"
            / f"ses-{ses:02d}"
            / task
        )

        icthresh = run_root / "melodic" / "ICthresh_zstat.nii.gz"
        timecourse = run_root / "melodic" / "report" / f"t{ic}.txt"

        if not icthresh.is_file():
            raise RuntimeError(f"Missing IC image source: {icthresh}")
        if not timecourse.is_file():
            raise RuntimeError(f"Missing IC time course: {timecourse}")

        out_map = blinded_dir / f"{blind_id}_ICmap.nii.gz"
        out_ts = blinded_dir / f"{blind_id}_timeseries.txt"

        subprocess.run(
            ["fslroi", str(icthresh), str(out_map), str(ic - 1), "1"],
            check=True,
        )
        shutil.copy2(timecourse, out_ts)

        scoring_rows.append({
            "blind_id": blind_id,
            "visual_label": "",
            "confidence": "",
            "notes": "",
        })

        key_rows.append({
            "blind_id": blind_id,
            "sub_id": sub,
            "ses_id": ses,
            "task_name": task,
            "PotentialICs": ic,
            "reference_label_source": row["gold_label_source"],
            "consensus_reference_label": int(float(row["manual_label"])),
            "historical_CICADA_label": int(float(row["baseline_auto_label"])),
            "candidate_11mm_CICADA_label": int(float(row["corrected_auto_label"])),
            "consensus_error_corrected": int(float(row["manual_error_corrected"])),
            "consensus_error_introduced": int(float(row["manual_error_introduced"])),
            "historical_smoothing_retention": row["historical_smoothing_retention"],
            "candidate_11mm_smoothing_retention": row["corrected_smoothing_retention"],
            "delta_smoothing_retention": row["delta_smoothing_retention"],
        })

    write_csv(
        blinded_dir / "BLINDED_SCORING.csv",
        scoring_rows,
        ["blind_id", "visual_label", "confidence", "notes"],
    )

    write_csv(
        key_dir / "PRIVATE_UNBLINDING_KEY.csv",
        key_rows,
        list(key_rows[0].keys()),
    )

    readme = f"""CICADA LOCKED HELD-OUT BLINDED IC REVIEW
========================================

PRIVATE HUMAN-SUBJECT-DERIVED MATERIAL.
DO NOT COMMIT OR UPLOAD TO PUBLIC GITHUB.

Changed ICs included: {len(flipped)}
Randomization seed: {seed}

Purpose
-------
These are ONLY the ICs whose final CICADA signal/noise classification differed
between:
  - historical CICADA smoothing: sigma=6 mm (~14.13 mm FWHM)
  - frozen candidate: 11 mm FWHM

The comparison was already performed on the locked KD/LS/MM consensus cohort.
This visual review is qualitative adjudication of the changed cases; it is NOT
parameter tuning.

Blinded review
--------------
Inspect each Review_### spatial map without opening the unblinding key.

Suggested visual_label:
  signal
  noise
  uncertain

Suggested confidence:
  high
  medium
  low

Optional:
  inspect Review_###_timeseries.txt if useful.

For consistent spatial viewing, use a standard MNI underlay such as:
  $FSLDIR/data/standard/MNI152_T1_2mm_brain.nii.gz

Important
---------
Complete BLINDED_SCORING.csv BEFORE opening:
  ../unblinding/PRIVATE_UNBLINDING_KEY.csv

Do not choose or test another smoothing kernel based on these judgments.
"""

    (blinded_dir / "README_BLINDED.txt").write_text(
        readme, encoding="utf-8"
    )

    print("=" * 76)
    print("CICADA HELD-OUT BLINDED IC REVIEW PACKET")
    print("=" * 76)
    print(f"Changed ICs included: {len(flipped)}")
    print(f"Randomization seed: {seed}")
    print(f"Blinded review folder:\n  {blinded_dir}")
    print(f"Private unblinding key:\n  {key_dir / 'PRIVATE_UNBLINDING_KEY.csv'}")
    print("\nDo NOT open the unblinding key until BLINDED_SCORING.csv is complete.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
