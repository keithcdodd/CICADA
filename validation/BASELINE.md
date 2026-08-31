# CICADA Development Validation Baseline

Frozen: 2026-08-31

## Purpose

This baseline freezes CICADA behavior before the 2026-2027 scientific
hardening and validation work. Candidate scientific changes should be
evaluated against this baseline using fixed validation datasets and
manual IC labels.

## Repository baseline

Repository: keithcdodd/CICADA
Branch at freeze: main

Commit:
8afc2e28350ef28d380b9ab286eac33a9c246073

Short commit:
8afc2e2

Git description:
v1.1.0-22-g8afc2e2

Validation baseline tag:
validation-baseline-2026-08-31

Historical released reference:
v1.1.0

## Software environment

FSL:
6.0.5.2
Build/commit: dc6f4207

MATLAB:
24.2.0.2740171
R2024b Update 1

Operating system:
macOS 15.3
Build 24D60

Relevant MATLAB toolboxes:
- Image Processing Toolbox 24.2
- Statistics and Machine Learning Toolbox 24.2

Additional CICADA MATLAB dependencies will be audited and added here
if required.

## v1.1.0 to development-baseline check

The development baseline contains 22 commits after v1.1.0, primarily
software robustness, QC, support-product, startup, comparison-file,
and wrapper consistency fixes.

The post-v1.1.0 modification to
templates/resampled_template.nii.gz was specifically checked.

Comparison with v1.1.0:
- NIfTI voxel data are identical.
- Voxelwise subtraction range: 0 to 0.
- Difference mean: 0.
- Difference standard deviation: 0.
- Nonzero voxel count and image summary statistics are identical.
- The only observed fslhd difference, aside from filename, is the
  NIfTI `descrip` metadata field reflecting the FSL version that wrote
  the file.

Therefore, the template modification does not represent a change in
image voxel values.

## Validation data

Original validation datasets:
TO BE POPULATED DURING ROADMAP ITEM 0B.

Manual gold-standard IC labels:
TO BE POPULATED DURING ROADMAP ITEM 0B.

Reference baseline outputs:
TO BE GENERATED DURING ROADMAP ITEM 0B.

## Baseline policy

Scientific behavior-changing experiments should compare against this
frozen baseline unless explicitly stated otherwise.

The historical v1.1.0 release remains available for reproduction of
the released CICADA implementation.
