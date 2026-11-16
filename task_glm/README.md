# CICADA Task GLM — optional standalone companion

**Version:** `0.3.3-experimental` (candidate for CICADA v1.4.0). This optional companion remains experimental; the historical v0.3.2-rc1 validation checkpoint below does not validate every later change.

**Entry points:** `cicada_task_glm.m` (preflight/fit) and `cicada_task_glm_export.m` (export existing fit).  
**Scope:** one completed single-run CICADA task scan, in two distinct modes: `original_joint` (default; original pre-ICA-cleanup BOLD with jointly fitted task, drift and selected ICA-noise regressors) or **experimental** `group_matched` (producer-verified, unsmoothed/unscaled Group-prepared BOLD with matched temporal operator and rank-aware covariance). Neither mode reruns ICA, CICADA denoising or Group QC.

**Outputs:** task COPE, **provisional** plug-in VARCOPE, AR covariance-parameter maps, valid/attempted masks, exact design, hashes and provenance. **No p/t/z or calibrated confidence intervals.**

## What users run

Add **only** the `task_glm/` folder to the MATLAB path; do not `genpath` or permanently add `helper_functions/`. The main entry points temporarily load the helper folder and restore the prior path. `tests/` is for developers, not routine analysis.

```matlab
addpath('/path/to/CICADA/task_glm','-begin');
which cicada_task_glm -all
which cicada_task_glm_export -all
```

Use a NEW analysis-output path in an already-existing parent outside CICADA-managed derivatives. Existing CICADA products are read-only; no Auto, Manual, Group, MELODIC, or FILM runs are triggered. The `cicada_task_glm_preflight.m` helper performs the read-only design check for **both** preflight-only and fitting calls; it does not itself execute FILM or fit BOLD. **A plain `cicada_task_glm(cfg)` call defaults to fitting and writing results. Always set `cfg.dryRun=true` for a preflight-only call.**

### Single-run named contrast — preflight only

```matlab
cfg=struct();
cfg.scanDir='/path/to/derivatives/cicada/sub-XXX/ses-01/foodpics_run-01';
cfg.outputDir='/path/to/separate_analysis_root/new_result_directory';
cfg.namedContrasts=struct('name','H_minus_N', ...
    'conditions',{{'highhedonic','nonfoods'}},'weights',[1,-1]);
cfg.icSelection='auto';    % 'manual' after completed Manual CICADA
cfg.driftModel='cosine'; cfg.highPassHz=0.008;
cfg.model='AR1';           % explicitly AR1 or AR2; no automatic selection
cfg.dryRun=true;          % REQUIRED for preflight only (no fit or writes)
R=cicada_task_glm(cfg);   % read-only design validation
```

The named condition labels above are **foodpics examples**, not hardcoded task requirements. The library consumes the scan's `hrf_conditions_raw.csv` labels. For more contrasts, provide a `cfg.namedContrasts` struct array. A genuinely different task/run is **not yet independently validated**; check events, task-HRF timing, label integrity, and design rank before using results.

### Fit without writing

```matlab
cfg.dryRun=false;          % reset after the preflight example above
cfg.runAnalysis=true; cfg.writeOutputs=false;
cfg.maxVoxels=256; cfg.batchVoxels=64;
R=cicada_task_glm(cfg);    % in-memory fit plus export preview; no files written
```

To explicitly fit the entire CICADA processing mask, use `cfg.maxVoxels=Inf; cfg.allowAllVoxels=true; cfg.allowLargeRun=true;`. A full 4D compressed BOLD image is loaded into memory before statistical batching. Validate resources first.

### Export a completed fit without refitting

```matlab
o=struct('writeOutputs',false); % preview only
P=cicada_task_glm_export(R,o);
% To actually write: NEW output directory only, after review:
o.writeOutputs=true;
% For sparse sampled fits ONLY: o.allowSparseTrial=true;
P=cicada_task_glm_export(R,o);
```

For a one-call fit and export to a new directory, explicitly set `cfg.dryRun=false; cfg.runAnalysis=true; cfg.writeOutputs=true` (or use the normal fit-and-export defaults in a fresh `cfg` without `dryRun` or legacy execution switches). Sparse samples need `allowSparseTrial=true` for export, use `sparseTrial_` filenames, and store NaN outside fitted voxels. Complete mask attempts retain attempted/valid distinctions; invalid voxels remain NaN. Each export is staged and re-read to verify geometry and values. No-overwrite protection is enforced.

### Experimental Group-matched mode — preflight before fitting

This mode is **not** the same as fitting arbitrary cleaned BOLD. It requires a
producer-time, hash-verified Group recipe from the **same** completed scan and
ICA classification, with `matched_temporal_ica`, degree-2 polynomial drift,
zero smoothing, and no voxelwise scaling. The current code does **not**
automatically support older Group output lacking matching producer attestation,
legacy sequential processing, spatial smoothing, or a different temporal recipe.
Do not retroactively label historical Group derivatives as verified.

```matlab
cfg=struct();
cfg.scanDir='/path/to/derivatives/cicada/sub-XXX/ses-01/foodpics_run-01';
cfg.outputDir='/path/to/separate_analysis_root/new_group_result';
cfg.namedContrasts=struct('name','H_minus_N', ...
    'conditions',{{'highhedonic','nonfoods'}},'weights',[1,-1]);
cfg.inputMode='group_matched';
cfg.groupRecipe='/path/to/producer_attested_group_recipe.json';
cfg.driftModel='polynomial'; cfg.driftDegree=2;
cfg.model='AR1';                 % or 'AR2'; no auto model selection
cfg.dryRun=true;
R=cicada_task_glm(cfg);       % verifies Group lineage and design, no writes
```

After a successful preflight, a bounded in-memory fit may set
`cfg.dryRun=false; cfg.runAnalysis=true; cfg.writeOutputs=false;`
and choose `cfg.maxVoxels` / `cfg.batchVoxels`. Group NIfTI export
requires the additional, **explicit** setting `cfg.allowGroupExport=true`
with `cfg.writeOutputs=true`, and only to a new unprotected destination.
Group maps include COPE and **provisional** plug-in VARCOPE, never calibrated
p/t/z or confidence interval maps. For the GLM alone,
`cfg.analysisMaskFile` may name a strictly aligned binary analysis subset;
this does not change the ICA/Group processing mask or temporal operator.

The Group-matched implementation was numerically verified for foodpics
`sub-135` using AR(1)/AR(2) and complete broad-mask NIfTI export;
other tasks/TRs and general inferential calibration remain unverified.

## Scientific input contract

Requires completed compatible CICADA Auto or Manual outputs: original BOLD `funcfile.nii.gz`, 3D mask, MELODIC mixing matrix, noise classification, `hrf_conditions_raw.csv`, cleaned provenance JSON with task timing and IC-selection metadata, and corresponding cleaned output. The original BOLD can include preprocessing and despiking that happened before ICA regression. Task HRFs must be generated and aligned to **that scan's** BOLD sampling reference; the GLM companion does **not** generate HRFs from event TSVs or assume nominal onset time. The default `original_joint` numerical design includes an intercept, raw task regressors, chosen drift, and selected noise-IC time series fitted jointly. `group_matched` instead verifies Group-prepared BOLD and the exact recorded temporal recipe, works in a reduced temporal observation space, and propagates each fitted covariance model through the matched temporal operator. Fixed AR(1) or AR(2) covariance can be requested; parameters are estimated by restricted likelihood on a finite grid. Run-specific design rank/estimability is checked. One contiguous run only; no multi-run boundaries or censoring.

### Dependencies

MATLAB with standard matrix routines, `niftiread`, `niftiinfo`, `niftiwrite`, table/JSON tools and Java available for canonical path, SHA-256 and staging. FSL is **not required** for the standalone REML analysis. Historical FILM comparisons used external tooling; the formerly bundled experimental FILM pilot was removed. The current `cicada_task_glm` path does not run FILM. Input files may be large; SHA-256 is read in chunks but full 4D `niftiread` is not streamed.

## Reproducibility and limits

New original-joint fits record SHA-256 hashes for at least **six input files** and **nine code files**, plus MATLAB and companion versions. An independent `analysisMaskFile` adds an input fingerprint; Group-matched mode also fingerprints its verified Group-prepared BOLD, producer recipe and additional preflight code. Before an authorized write, current hashes are compared with the pre-fit fingerprints; no hashing is required to generate preflight-only designs. The version stored with a new fit must match the exact companion source used to generate it; new fits from this candidate record `0.3.3-experimental`; verify the exact source hashes in each run's provenance. Existing v0.3.1 and v0.3.2 outputs retain their original version provenance and must not be relabeled. Stronger provenance does not validate the scientific correctness of inputs.

- **No inferential significance output.** Plug-in VARCOPE does not incorporate full uncertainty due to covariance-parameter estimation or adaptive ICA estimation and task-informed classification. Simulation experiments found coverage deficits and false-positive inflation under some covariance misspecification/nonstationary conditions. Do not treat VARCOPE as calibrated for p-values, t/z maps or confidence intervals.
- The processing mask is not an analysis-quality mask. Additional mean-intensity, tSNR and anatomical/geometric QC are important, especially near DBS electrodes.
- Manual preflight and bounded fitting were tested, but the compared Auto and Manual runs selected the **same 52 noise IC indices**. Behavior when Manual classifications actually differ has not been independently exercised.
- Whole-mask fitting and export were validated using Auto CICADA on two foodpics subjects (sub-132 and sub-135). The local inventory of raw-HRF files contained three scan directories, all `foodpics_run-01`; a different task/TR and full end-to-end adaptive-ICA statistical calibration remain untested.
- Output uses single-precision spatial maps; analysis design and model settings are preserved separately in `analysis_model.mat` and `provenance.json`.

## Validation record and release criteria

See `VALIDATION_SUMMARY.md` for the historical original-joint numerical and MATLAB tests and `RELEASE_CHECKLIST.md` for remaining external-publication checks. Later Group-matched tests are summarized above and in the integrated CICADA development handoff. The v0.3.2-rc1 bounded MATLAB smoke test passed before these later development changes. This package is **eligible for an explicitly experimental companion release after owner review of the final Git diff**, but is not yet officially released. An altered-Manual-label fixture, different task/TR, final clean-install compatibility, and end-to-end adaptive ICA inference validation remain future scope; these are **not** claimed as tested.

## Licensing

The parent CICADA repository carries the **GNU GPL v3**. This optional companion is distributed with a copy of the same `LICENSE` for independent distribution. Preserve license notices, attribution and source distribution requirements as appropriate. See the parent repository for any other applicable notices.
