# Option 1 standalone companion — release gate (v0.3.2-rc1)

**Classification:** candidate for an **optional experimental, noninferential companion**, not a fully validated general task fMRI inference engine. No production CICADA processing modifications.

## Passed on user's MATLAB

- [x] AR1/AR2 numerical reference tests and direct GLS/REML objective/optimizer checks.
- [x] Two-subject Auto whole-brain AR1 fits, mask/invalid voxel accounting and spatial export read-back.
- [x] Manual preflight and bounded model fit.
- [x] Three simultaneous contrasts via standalone entry point; separate export preview; path restored.
- [x] Input lineage / output overwrite / managed-tree protection checks.
- [x] SHA256 known fixture, pre-fit six input + nine source fingerprints, 16-voxel provenance export with 11 verified files (v0.3.1).
- [x] Fixed-ICA null, effect recovery, covariance misspecification and nonstationarity validation performed; limitations recorded, NOT resolved for inference.

## Required before an experimental companion release

- [x] MATLAB v0.3.2-rc1 smoke test on sub-135: 4/4 valid, expected version string, `PREVIEW_ONLY`, no files written.
- [x] Bundled `LICENSE` was compared byte-for-byte against `CICADA/LICENSE` in the submitted repository ZIP (identical GPLv3 text); user guide and example reviewed statically. The exact minimum MATLAB release is not asserted because it has not been independently established.
- [ ] Before external publication, confirm any contributor attributions/notices and test the example on a clean MATLAB path if claiming clean-install compatibility.
- [x] Scope defined: **foodpics run-01 with matching raw HRF and user-verified task timing**, original-BOLD single-run Auto ICA model and bounded Manual branch, exploratory noninferential output. Other tasks/TRs are not validated.
- [x] Auto and Manual noise-IC index sets compared: both contain the **same 52 indices**; branch compatibility, **not** changed-label behavior, is supported by these tests.
- [ ] Optional future validation: change noise-IC selection in a separate synthetic/sandbox fixture if a stronger Manual-difference claim is desired; never change managed scan classifications for testing.
- [ ] Confirm status of existing repository wrapper edits before commit; commit only `task_glm/` for the companion when ready (or intentionally separate other staged changes), no push/release without user's authorization.

## Experimental-release decision (software gate)

The v0.3.2-rc1 MATLAB smoke and the documentation/scope review are complete. **Suitable for a user-authorized local commit or narrowly labeled public experimental companion release**, conditional on checking repository status/contributor notices and reviewing the final Git diff. This is not an inferential validation or evidence of compatibility across arbitrary tasks.

## Separate research track for inferential claims

- [ ] Calibrate adaptive MELODIC extraction, task-informed ICA classification, and parameter selection under null and injected signals, using isolated synthetic CICADA processing trees.
- [ ] Assess realistic nuisance motion, physiological artifacts, temporal nonstationarity, and spatial multiple-testing behavior.
- [ ] Validate a prespecified inferential uncertainty/calibration method; no p/t/z or CI maps before this gate.
- [ ] Different task/TR and run-boundary/censoring features require new science and software tests before claiming support.

## Release decision

This package remains **unreleased** until the repository owner chooses to commit, push, or tag it. A user-approved experimental companion release may state the validated foodpics scope, with **no inferential significance output** and explicit limitations. General task support and fully calibrated inference are separate future work.
