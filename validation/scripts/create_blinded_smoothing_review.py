#!/usr/bin/env python3
"""
Create a PRIVATE, BLINDED visual-review packet for CICADA ICs whose final
automatic classification changed in the FIX-training smoothing experiment.

Inputs
------
--comparison:
    fixtraining_smoothing_ic_level_PRIVATE.csv
--data-root:
    CICADA derivatives root containing sub-* folders
--output-dir:
    PRIVATE directory OUTSIDE Git

Outputs
-------
<output-dir>/blinded_review/
    Review_001_ICmap.nii.gz
    Review_001_timeseries.txt
    ...
    BLINDED_SCORING.csv
    README_BLINDED.txt

<output-dir>/unblinding/
    PRIVATE_UNBLINDING_KEY.csv

The blinded review folder contains neutral IDs only. The unblinding key stores
subject/task/IC identity, historical KD label, baseline CICADA label, corrected
label, and whether the original quantitative comparison classified the flip
as an error corrected or introduced.

All outputs remain human-subject-derived and PRIVATE.
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path
import random
import secrets
import shutil
import subprocess
import sys


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--comparison", required=True, type=Path)
    p.add_argument("--data-root", required=True, type=Path)
    p.add_argument("--output-dir", required=True, type=Path)
    p.add_argument(
        "--seed",
        type=int,
        default=None,
        help="Optional randomization seed. Default: securely generated.",
    )
    p.add_argument("--allow-output-inside-git", action="store_true")
    return p.parse_args()


def find_git_root(path: Path):
    path = path.resolve()
    for candidate in (path, *path.parents):
        if (candidate / ".git").exists():
            return candidate
    return None


def require_program(name):
    path = shutil.which(name)
    if path is None:
        raise RuntimeError(f"Required program not found on PATH: {name}")
    return path


def write_csv(path, rows, fields):
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def main():
    args = parse_args()

    args.comparison = args.comparison.expanduser().resolve()
    args.data_root = args.data_root.expanduser().resolve()
    args.output_dir = args.output_dir.expanduser().resolve()

    git_root = find_git_root(args.output_dir)
    if git_root is not None and not args.allow_output_inside_git:
        raise RuntimeError(
            "\nPRIVACY SAFETY STOP\n"
            "Participant-derived blinded review output was requested inside "
            "a Git repository.\n\n"
            f"Requested output: {args.output_dir}\n"
            f"Git repository:    {git_root}\n"
        )

    require_program("fslroi")

    if not args.comparison.is_file():
        raise RuntimeError(f"Comparison CSV not found: {args.comparison}")

    with args.comparison.open(newline="", encoding="utf-8-sig") as f:
        reader = csv.DictReader(f)
        rows = list(reader)
        fields = set(reader.fieldnames or [])

    required = {
        "sub_id",
        "ses_id",
        "task_name",
        "PotentialICs",
        "manual_label",
        "baseline_auto_label",
        "corrected_auto_label",
        "label_flipped",
        "manual_error_corrected",
        "manual_error_introduced",
        "historical_smoothing_retention",
        "corrected_smoothing_retention",
        "delta_smoothing_retention",
    }
    missing = required - fields
    if missing:
        raise RuntimeError(
            "Comparison CSV missing required columns: "
            + ", ".join(sorted(missing))
        )

    flipped = [
        r for r in rows
        if int(float(r["label_flipped"])) == 1
    ]

    if not flipped:
        raise RuntimeError("No final-label-flipped ICs found.")

    seed = args.seed if args.seed is not None else secrets.randbits(63)
    rng = random.Random(seed)
    rng.shuffle(flipped)

    blinded_dir = args.output_dir / "blinded_review"
    key_dir = args.output_dir / "unblinding"

    if blinded_dir.exists() or key_dir.exists():
        raise RuntimeError(
            "Output review directories already exist. "
            "Choose a new --output-dir so a previous blinded review is not "
            "silently overwritten."
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
            [
                "fslroi",
                str(icthresh),
                str(out_map),
                str(ic - 1),
                "1",
            ],
            check=True,
        )

        shutil.copy2(timecourse, out_ts)

        scoring_rows.append(
            {
                "blind_id": blind_id,
                "visual_label": "",
                "confidence": "",
                "notes": "",
            }
        )

        key_rows.append(
            {
                "blind_id": blind_id,
                "sub_id": sub,
                "ses_id": ses,
                "task_name": task,
                "PotentialICs": ic,
                "historical_KD_manual_label": int(float(row["manual_label"])),
                "baseline_CICADA_label": int(float(row["baseline_auto_label"])),
                "corrected_CICADA_label": int(float(row["corrected_auto_label"])),
                "manual_error_corrected":
                    int(float(row["manual_error_corrected"])),
                "manual_error_introduced":
                    int(float(row["manual_error_introduced"])),
                "historical_smoothing_retention":
                    row["historical_smoothing_retention"],
                "corrected_smoothing_retention":
                    row["corrected_smoothing_retention"],
                "delta_smoothing_retention":
                    row["delta_smoothing_retention"],
            }
        )

    write_csv(
        blinded_dir / "BLINDED_SCORING.csv",
        scoring_rows,
        ["blind_id", "visual_label", "confidence", "notes"],
    )

    write_csv(
        key_dir / "PRIVATE_UNBLINDING_KEY.csv",
        key_rows,
        [
            "blind_id",
            "sub_id",
            "ses_id",
            "task_name",
            "PotentialICs",
            "historical_KD_manual_label",
            "baseline_CICADA_label",
            "corrected_CICADA_label",
            "manual_error_corrected",
            "manual_error_introduced",
            "historical_smoothing_retention",
            "corrected_smoothing_retention",
            "delta_smoothing_retention",
        ],
    )

    readme = f"""CICADA blinded visual review
============================

PRIVATE HUMAN-SUBJECT-DERIVED MATERIAL.
DO NOT COMMIT OR UPLOAD TO PUBLIC GITHUB.

Number of changed ICs: {len(flipped)}

Randomization seed:
{seed}

Review instructions
-------------------
For each neutral Review_### case, inspect the spatial IC map and associated
time course WITHOUT opening the unblinding key.

Suggested visual_label values:
    signal
    noise
    uncertain

Suggested confidence values:
    high
    medium
    low

Do not use the historical KD label, baseline CICADA label, corrected CICADA
label, subject ID, task identity, or whether the flip was previously counted
as "corrected" versus "introduced" until all visual classifications are saved.

Spatial maps
------------
Review_###_ICmap.nii.gz is the corresponding 3-D volume extracted from
ICthresh_zstat.nii.gz.

For consistent viewing, load a standard MNI underlay in FSLeyes, for example:
    $FSLDIR/data/standard/MNI152_T1_2mm_brain.nii.gz

Time courses
------------
Review_###_timeseries.txt contains the corresponding MELODIC IC time course.

After review
------------
Complete BLINDED_SCORING.csv first. Then the separate private key can be used
to compare the new blinded visual classification against:
    - historical KD manual classification
    - baseline CICADA
    - corrected 6-mm-FWHM CICADA

Randomization metadata belongs only in this PRIVATE review package.
"""

    (blinded_dir / "README_BLINDED.txt").write_text(
        readme, encoding="utf-8"
    )

    print("=" * 72)
    print("CICADA BLINDED IC REVIEW PACKET")
    print("=" * 72)
    print(f"Changed ICs included: {len(flipped)}")
    print(f"Randomization seed: {seed}")
    print(f"Blinded review folder:\n  {blinded_dir}")
    print(f"Private unblinding key:\n  {key_dir / 'PRIVATE_UNBLINDING_KEY.csv'}")
    print("\nDo NOT open the unblinding key until BLINDED_SCORING.csv is complete.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
