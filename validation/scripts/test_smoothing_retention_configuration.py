#!/usr/bin/env python3
"""
Static regression checks for CICADA's Smoothing_Retention kernel configuration.

This test uses only source code text; it contains no participant-derived data.
Run from the CICADA repository root:

    python validation/scripts/test_smoothing_retention_configuration.py
"""

from pathlib import Path
import math
import re

SOURCE = Path("basescripts/CICADA_1_MasksandICAs.sh")
FWHM_PER_SIGMA = 2.354820045


def main():
    if not SOURCE.exists():
        raise SystemExit(f"Cannot find {SOURCE}; run from repository root.")

    text = SOURCE.read_text(encoding="utf-8")

    assert 'CICADA_SMOOTHING_RETENTION_MODE="${CICADA_SMOOTHING_RETENTION_MODE:-revised}"' in text
    assert 'CICADA_SMOOTHING_RETENTION_FWHM_MM="11"' in text
    assert 'CICADA_SMOOTHING_RETENTION_SIGMA_MM="4.671269902"' in text
    assert 'CICADA_SMOOTHING_RETENTION_FWHM_MM="14.12892027"' in text
    assert 'CICADA_SMOOTHING_RETENTION_SIGMA_MM="6"' in text

    expected_sigma = 11.0 / FWHM_PER_SIGMA
    assert math.isclose(expected_sigma, 4.671269902, rel_tol=0, abs_tol=5e-10)

    revised_cmd = (
        'fslmaths ${mel_fol}/ICthresh_zstat.nii.gz '
        '-s "${CICADA_SMOOTHING_RETENTION_SIGMA_MM}" -abs '
        '"${ROIcalcfol}/fullvolICA_tmp_smoothed_nothresh.nii.gz"'
    )
    assert text.count(revised_cmd) == 1, (
        "Expected exactly one parameterized Smoothing_Retention fslmaths command."
    )

    # The old hard-coded Smoothing_Retention command must be gone.
    old_exact = (
        'fslmaths ${mel_fol}/ICthresh_zstat.nii.gz -s 6 -abs '
        '"${ROIcalcfol}/fullvolICA_tmp_smoothed_nothresh.nii.gz"'
    )
    assert old_exact not in text

    # Other -s 6 operations are allowed and intentionally not globally changed.
    other_sigma6 = len(re.findall(r'(?<![0-9.])-s\s+6(?![0-9.])', text))
    print("PASS: revised 11-mm-FWHM Smoothing_Retention configuration present.")
    print("PASS: historical sigma=6 reproducibility mode present.")
    print("PASS: old hard-coded Smoothing_Retention `-s 6` command removed.")
    print(f"INFO: {other_sigma6} other `-s 6` occurrence(s) remain in this source by design.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
