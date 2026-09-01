#!/usr/bin/env python3
"""
Aggregate unblinding for the CICADA locked held-out 11-mm smoothing review.

Reads the completed blinded scoring CSV and the separate private unblinding key.
Prints aggregate agreement only; does not print subject IDs, task names, or IC
numbers.

All material remains PRIVATE human-subject-derived validation data.
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--scoring", required=True, type=Path)
    p.add_argument("--key", required=True, type=Path)
    p.add_argument("--output", type=Path, default=None)
    return p.parse_args()


def read_by_id(path):
    with path.open(newline="", encoding="utf-8-sig") as f:
        rows = list(csv.DictReader(f))
    if not rows:
        raise RuntimeError(f"No rows in {path}")
    out = {}
    for row in rows:
        bid = row["blind_id"].strip()
        if bid in out:
            raise RuntimeError(f"Duplicate blind_id: {bid}")
        out[bid] = row
    return out


def visual_value(v):
    x = v.strip().lower()
    if x == "signal":
        return 1
    if x == "noise":
        return 0
    if x in {"uncertain", "unsure", "?"}:
        return None
    raise RuntimeError(f"Unexpected visual_label: {v!r}")


def main():
    args = parse_args()
    scoring = read_by_id(args.scoring.expanduser().resolve())
    key = read_by_id(args.key.expanduser().resolve())

    if set(scoring) != set(key):
        raise RuntimeError("Scoring and key blind-ID sets differ.")

    rows = []
    for bid in sorted(scoring):
        s = scoring[bid]
        k = key[bid]
        rows.append({
            "blind_id": bid,
            "visual": visual_value(s["visual_label"]),
            "confidence": s.get("confidence", "").strip().lower(),
            "consensus": int(float(k["consensus_reference_label"])),
            "historical": int(float(k["historical_CICADA_label"])),
            "candidate": int(float(k["candidate_11mm_CICADA_label"])),
        })

    certain = [r for r in rows if r["visual"] is not None]
    uncertain = [r for r in rows if r["visual"] is None]

    for r in rows:
        if r["historical"] == r["candidate"]:
            raise RuntimeError("Key contains a case that was not a final-label flip.")

    n = len(certain)
    sig = sum(r["visual"] == 1 for r in certain)
    noise = sum(r["visual"] == 0 for r in certain)

    hist_agree = sum(r["visual"] == r["historical"] for r in certain)
    cand_agree = sum(r["visual"] == r["candidate"] for r in certain)
    cons_agree = sum(r["visual"] == r["consensus"] for r in certain)

    lines = [
        "=" * 76,
        "CICADA HELD-OUT 11-MM BLINDED REVIEW — AGGREGATE UNBLINDING",
        "=" * 76,
        f"Reviewed changed ICs: {len(rows)}",
        f"Definite visual labels: {n}",
        f"Uncertain: {len(uncertain)}",
        "",
        "Fresh blinded visual classifications",
        f"  Signal: {sig}",
        f"  Noise:  {noise}",
        "",
        "Agreement with fresh blinded visual review",
        f"  Historical CICADA: {hist_agree}/{n}",
        f"  11-mm CICADA:      {cand_agree}/{n}",
        f"  Consensus label:   {cons_agree}/{n}",
        "",
        "Among the smoothing-induced final-label flips",
        f"  Fresh visual review favors historical: {hist_agree}",
        f"  Fresh visual review favors 11 mm:      {cand_agree}",
        f"  Uncertain:                           {len(uncertain)}",
        "",
        "By confidence",
    ]

    # Preserve confidence labels exactly as entered, including "very low".
    seen = []
    for r in certain:
        c = r["confidence"] or "unspecified"
        if c not in seen:
            seen.append(c)

    for c in seen:
        subset = [r for r in certain if (r["confidence"] or "unspecified") == c]
        h = sum(r["visual"] == r["historical"] for r in subset)
        cand = sum(r["visual"] == r["candidate"] for r in subset)
        cons = sum(r["visual"] == r["consensus"] for r in subset)
        lines.append(
            f"  {c}: n={len(subset)}, historical={h}, "
            f"11mm={cand}, consensus={cons}"
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
    raise SystemExit(main())
