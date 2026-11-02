# CICADA v1.3.0

**Published:** October 26, 2026  
**Release tag:** `v1.3.0` (`3aba2344a03bd6f6032ba3164cd1b182621446c3`)  
**Comparison:** [`v1.2.0...v1.3.0`](https://github.com/keithcdodd/CICADA/compare/v1.2.0...v1.3.0)

CICADA v1.3.0 improves **processing provenance, task-design timing, and reproducibility** for Automatic, Manual, and Group CICADA. This release does not add a task-GLM analysis pipeline; task-GLM work developed after v1.3.0 is separate.

## Added and improved

### Processing provenance

- Automatic and Manual CICADA write processing-provenance JSON sidecars for their cleaned BOLD outputs, recording source inputs, IC-selection information, denoising choices, available task-HRF information, and source Git revision metadata.
- Group CICADA records processing decisions and scan/source information in a standalone `group_qc_provenance.json` and its MATLAB QC structure, including relevant temporal and spatial processing details.
- The provenance is intended to make the actual processing configuration auditable; it does not itself establish statistical validity.

### Task timing and HRF regressors

- Added `cicada_resolve_task_timing.m` to resolve task-HRF sampling times from an explicit user/pipeline reference, an appropriate **preprocessed-BOLD JSON** sidecar, or a documented fallback when the reference cannot be verified.
- Added `cicada_task_hrf_design.m` for continuous-time, event-based task-HRF design construction at the resolved acquisition/reference sample times.
- Automatic CICADA now saves **unstandardized condition regressors** as `hrf_conditions_raw.csv` for downstream analysis, while retaining the normalized `hrf_conditions.csv` (including its derived `Combined` column) for existing IC-classification/QC purposes.
- Added fractional timing and timing-resolution validation scripts. Task timing still requires scrutiny against each scan's preprocessing history; unverified fallbacks should not be interpreted as confirmed acquisition timing.

### File handling and compatibility

- Stage 1 preserves the matching fMRIPrep preprocessed-BOLD JSON as `funcfile_preproc_bold.json` when available and avoids retaining stale metadata when unavailable.
- Added file-discovery logic that ignores hidden macOS AppleDouble (`._*`) files, reducing ambiguous matches and accidental selection of metadata stubs.
- Added an option in Stage 1 for a user-provided subependymal anatomical probability map, in place of the generated default when explicitly supplied.

### Reproducible quality-control sampling

- Subject-level and group-level random samples used in QC have deterministic scan-related seeds while restoring MATLAB's prior RNG state afterward.
- Corrected group-level IC ranking logic to identify the lowest-scoring k-means cluster by its **centroid**, rather than assuming a stable numeric cluster label.

## Continued behavior and scope

- Automatic and Manual CICADA remain the principal subject-level ICA classification/denoising workflows; Group CICADA remains a separate group processing/QC stage.
- Group CICADA's matched temporal cleanup default and optional bandpass behavior were introduced in **v1.2.0**, not v1.3.0; see `RELEASE_NOTES_v1.2.0.md`.
- A separate `task_glm/` experimental companion was **not included in the v1.3.0 tag**. If released, it should be documented under its own later release or companion version.

## Validation and limitations

The repository adds `validation/scripts/test_cicada_resolve_task_timing.m` and `validation/scripts/test_cicada_hrf_fractional.m`. These targeted checks support the new task-timing and HRF behavior. They do not constitute universal validation of task-event alignment, alternative preprocessing configurations, or downstream inferential GLMs. Exact release behavior should be reproduced from the tagged source and recorded input provenance.

## Source comparison

The v1.3.0 release tag includes **five commits after v1.2.0** and changes to 18 tracked paths (including the newly committed v1.2.0 release-notes file). The relevant code additions span the Automatic/Manual/Group pipelines and helpers. These notes describe changes present at the **v1.3.0 tag**; subsequent main-branch commits are not part of this release.

See also: [CICADA repository](https://github.com/keithcdodd/CICADA) and `CICADA_UserGuide.docx`.
