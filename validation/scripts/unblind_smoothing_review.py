#!/usr/bin/env python3
"""
Aggregate unblinding analysis for the CICADA smoothing blinded-review packet.

Reads:
  * BLINDED_SCORING.csv
  * PRIVATE_UNBLINDING_KEY.csv

Prints ONLY aggregate results by default. It does not print subject IDs,
task names, or IC numbers.

All inputs/outputs remain PRIVATE human-subject-derived validation material.
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path
import math
import sys


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--scoring", required=True, type=Path)
    p.add_argument("--key", required=True, type=Path)
    p.add_argument("--output", type=Path, default=None)
    return p.parse_args()


def read_by_id(path, id_col="blind_id"):
    with path.open(newline="", encoding="utf-8-sig") as f:
        rows = list(csv.DictReader(f))
    if not rows:
        raise RuntimeError(f"No rows in {path}")
    out = {}
    for row in rows:
        key = row[id_col].strip()
        if key in out:
            raise RuntimeError(f"Duplicate {id_col}: {key}")
        out[key] = row
    return out


def norm_visual(value):
    v = value.strip().lower()
    aliases = {
        "signal": 1,
        "noise": 0,
        "uncertain": None,
        "unsure": None,
        "?": None,
    }
    if v not in aliases:
        raise RuntimeError(
            f"Unexpected visual_label {value!r}; use signal/noise/uncertain"
        )
    return aliases[v]


def pct(num, den):
    return 100.0 * num / den if den else math.nan


def main():
    args = parse_args()
    scores = read_by_id(args.scoring)
    key = read_by_id(args.key)

    if set(scores) != set(key):
        raise RuntimeError(
            "Blind-ID sets differ between scoring file and unblinding key."
        )

    reviewed = []
    for blind_id in sorted(scores):
        s = scores[blind_id]
        k = key[blind_id]

        visual = norm_visual(s["visual_label"])
        confidence = s.get("confidence", "").strip().lower()

        reviewed.append({
            "blind_id": blind_id,
            "visual": visual,
            "confidence": confidence,
            "historical": int(float(k["historical_KD_manual_label"])),
            "baseline": int(float(k["baseline_CICADA_label"])),
            "corrected": int(float(k["corrected_CICADA_label"])),
        })

    n = len(reviewed)
    certain = [r for r in reviewed if r["visual"] is not None]
    uncertain = [r for r in reviewed if r["visual"] is None]

    visual_signal = sum(r["visual"] == 1 for r in certain)
    visual_noise = sum(r["visual"] == 0 for r in certain)

    baseline_agree = sum(r["visual"] == r["baseline"] for r in certain)
    corrected_agree = sum(r["visual"] == r["corrected"] for r in certain)
    historical_agree = sum(r["visual"] == r["historical"] for r in certain)

    # Because every case in this packet is a final-label flip, baseline and
    # corrected labels should differ. Verify this rather than assume it.
    nonflips = [
        r for r in reviewed if r["baseline"] == r["corrected"]
    ]
    if nonflips:
        raise RuntimeError(
            "Unblinding key contains cases where baseline and corrected "
            "labels do not differ."
        )

    baseline_favored = baseline_agree
    corrected_favored = corrected_agree

    by_conf = {}
    for conf in ("high", "medium", "low", ""):
        subset = [
            r for r in certain if r["confidence"] == conf
        ]
        if not subset:
            continue
        by_conf[conf or "unspecified"] = {
            "n": len(subset),
            "baseline": sum(
                r["visual"] == r["baseline"] for r in subset
            ),
            "corrected": sum(
                r["visual"] == r["corrected"] for r in subset
            ),
            "historical": sum(
                r["visual"] == r["historical"] for r in subset
            ),
        }

    lines = []
    lines.append("=" * 72)
    lines.append("CICADA BLINDED SMOOTHING REVIEW — AGGREGATE UNBLINDING")
    lines.append("=" * 72)
    lines.append(f"Reviewed changed ICs: {n}")
    lines.append(f"Definite visual labels: {len(certain)}")
    lines.append(f"Uncertain: {len(uncertain)}")
    lines.append("")
    lines.append("Fresh blinded visual classifications")
    lines.append(f"  Signal: {visual_signal}")
    lines.append(f"  Noise:  {visual_noise}")
    lines.append("")
    lines.append("Agreement with fresh blinded visual review")
    lines.append(
        f"  Baseline CICADA:  {baseline_agree}/{len(certain)} "
        f"({pct(baseline_agree, len(certain)):.1f}%)"
    )
    lines.append(
        f"  Corrected CICADA: {corrected_agree}/{len(certain)} "
        f"({pct(corrected_agree, len(certain)):.1f}%)"
    )
    lines.append(
        f"  Historical KD:    {historical_agree}/{len(certain)} "
        f"({pct(historical_agree, len(certain)):.1f}%)"
    )
    lines.append("")
    lines.append("Among the smoothing-induced final-label flips")
    lines.append(f"  Fresh visual review favors baseline:  {baseline_favored}")
    lines.append(f"  Fresh visual review favors corrected: {corrected_favored}")
    lines.append(f"  Uncertain:                         {len(uncertain)}")

    if by_conf:
        lines.append("")
        lines.append("By visual confidence")
        for conf, d in by_conf.items():
            lines.append(
                f"  {conf:11s}: n={d['n']}, "
                f"baseline={d['baseline']}, "
                f"corrected={d['corrected']}, "
                f"historical_KD={d['historical']}"
            )

    text = "\n".join(lines)
    print(text)

    if args.output is not None:
        out = args.output.expanduser().resolve()
        out.parent.mkdir(parents=True, exist_ok=True)
        out.write_text(text + "\n", encoding="utf-8")
        print(f"\nPrivate aggregate summary written to:\n  {out}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
