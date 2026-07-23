# CICADA v1.1.0

CICADA v1.1.0 adds additive subject-level support and reliability derivatives and organized group-level stacking while preserving legacy CICADA denoising, IC classification, Group CICADA outputs, and the legacy group signal mask.

## Added

- Direct functional-data reliability products.
- Robust tSNR and normalized functional-data QC products.
- Continuous retained-signal and retained-noise IC evidence maps.
- p50 and p95 thresholded signal/noise evidence derivatives.
- Smoothed and extended signal-evidence derivatives.
- Subject-level provenance and validation summaries under each task's `support_products/` folder.
- Group-level 4-D stacks and summary maps under:
  - `support_products/core/`
  - `support_products/qc/`
  - `support_products/provenance/`
- An authoritative group manifest containing zero-based and one-based volume indices, subject/session/task identifiers, CICADA mode, release metadata, task directories, and exact source paths.
- Staged group support-output generation, stack-volume validation, and atomic directory replacement.

## Primary downstream products

- `support_products/core/signal_ic_overlap_raw_maps.nii.gz`
- `support_products/core/direct_reliability_masks.nii.gz`
- `support_products/provenance/support_products_manifest.tsv`

Signal evidence is not a calibrated probability of neuronal signal. Thresholded and smoothed maps are QC and sensitivity derivatives, not definitive tissue masks.

## Compatibility

- No change to automatic or manual IC classification.
- No change to CICADA denoising calculations.
- No change to legacy subject outputs.
- No change to legacy Group CICADA output locations.
- No change to `group_funcmask.nii.gz` or `group_signal_funcmask.nii.gz`.

## Validation

The release candidate was validated with:

- multiple fresh automatic subject runs;
- one manual-path run;
- subject-level smoke tests and visual inspection;
- voxelwise legacy-output comparisons;
- a three-subject Group CICADA regression test;
- exact 4-D volume-order verification;
- voxel-identical comparisons between pre- and post-packaging support products;
- confirmation that legacy group masks remained unchanged.

## Intended downstream use

The continuous group signal-evidence and direct-reliability stacks provide the principal bridge to support-aware downstream methods, including IFFU. Analysis-specific eligibility masks and thresholds should be derived and validated within the downstream method rather than treated as a universal CICADA signal mask.
