# CICADA Smoothing_Retention kernel validation

## Decision

CICADA's revised default `Smoothing_Retention` kernel is **11 mm FWHM**
(Gaussian sigma = 4.671269902 mm for FSL `fslmaths -s`).

The historical CICADA behavior remains available as an explicit reproducibility
mode: sigma = 6 mm (FWHM approximately 14.12892 mm).

This smoothing operation is used to calculate the `Smoothing_Retention`
classifier feature. It is **not** conventional fMRI preprocessing smoothing and
is distinct from other CICADA smoothing operations used for masks or support
products.

## Why this was revisited

The original implementation used:

```bash
fslmaths ICthresh_zstat.nii.gz -s 6 -abs ...
```

FSL interprets `-s` as Gaussian sigma in millimeters. Thus the historical
operation used sigma = 6 mm, equivalent to approximately 14.13 mm FWHM,
rather than 6 mm FWHM.

Rather than directly replacing the historical operation with 6-mm FWHM, the
spatial scale of the feature was empirically reevaluated.

## Development-set calibration

The original FIX-training development scans (10 resting-state and 10 task
scans) were used for smoothing-scale development. Historical CICADA output was
first reconstructed exactly in a language-neutral Python regression harness.

Candidate smoothing scales demonstrated a bounded stable operating region:

- 6 mm FWHM: 10 final classification changes versus historical CICADA; 1 moved
  toward the historical KD manual label and 9 moved away.
- 8 mm FWHM: 5 classification changes; all 5 moved away from the KD manual
  label.
- 9, 10, 11, and 12 mm FWHM: no final classification changes.
- Historical sigma=6 mm / approximately 14.13 mm FWHM: reference behavior.
- 16 mm FWHM: deterioration began.
- 18-24 mm FWHM: progressively more adverse classification changes.

Isolated analysis of `Smoothing_Retention` also supported a central operating
region rather than either boundary. Eleven millimeters FWHM was frozen as a
central candidate before the held-out consensus cohort was examined.

A blinded review of the 10 development-set ICs changed by 6-mm FWHM favored
historical CICADA in 8/10 cases, supporting the conclusion that 6-mm FWHM was
too weak for the existing classifier.

## Locked held-out validation

The frozen 11-mm-FWHM candidate was then evaluated once on the original
consensus-labeled held-out schizophrenia cohort: 30 resting-state scans and
30 task scans, totaling 4,783 ICs.

Historical CICADA labels were reproduced exactly for all 60 scans before the
candidate comparison was accepted.

Aggregate performance was essentially preserved:

| Metric | Historical | 11 mm FWHM |
| --- | ---: | ---: |
| Mean scan-level accuracy | 0.9791 | 0.9790 |
| Mean scan-level signal F1 | 0.9235 | 0.9235 |
| Mean scan-level signal sensitivity | 0.9298 | 0.9310 |
| Mean scan-level noise sensitivity | 0.9922 | 0.9917 |

Only 4 of 4,783 final IC classifications changed. Relative to the historical
KD/LS/MM consensus labels, 2 changes corrected historical CICADA disagreements
and 2 introduced disagreements.

The four changed ICs were subsequently reassessed while blinded to the
historical CICADA classification, 11-mm classification, and consensus label.
Fresh visual/temporal review favored the 11-mm classification for all four
cases. This four-case review is supportive qualitative evidence and should not
be interpreted as an independent multi-rater validation study.

## Interpretation

The evidence does **not** establish that 11 mm FWHM is statistically superior
to the historical approximately 14.13-mm-FWHM behavior. Instead, it supports
11 mm FWHM as a deliberate, interpretable recalibration within an empirically
stable operating range that preserves held-out classification performance.

The historical implementation remains explicitly available for reproduction of
prior CICADA analyses.

## Privacy

Human-subject scan-level and IC-level validation outputs are private and are not
stored in the public repository. Public repository materials should contain only
code, documentation, synthetic/public fixtures, checksums where appropriate,
and explicitly approved aggregate summaries.
