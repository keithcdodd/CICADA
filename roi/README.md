# CICADA optional ROI time-series companion (experimental)

Entry point: `roi/cicada_roi_timeseries.m`. Add only the `roi/` folder to the MATLAB path, not `roi/helper_functions/` or `task_glm/helper_functions/`. The public function temporarily resolves its implementation and SHA-256 dependency and restores the original MATLAB path on return, including error return. It does not run or alter ICA, Group CICADA, task GLM, or MRI derivatives.

Input: one 4D BOLD NIfTI and any combination of explicitly typed 3D integer-labeled atlases, 3D/4D binary masks, or 3D/4D finite nonnegative weight maps. An ROI stack's fourth axis indexes ROIs, not time. Binary maps use ordinary means; weighted maps use nonnegative normalized weights after voxel QC. Fractional maps require explicit `roiTypes='weighted'`. Overlap across maps is allowed. Values need not sum to one before selection.

An aligned `cfg.qualityMaskFile`, or `cfg.scanDir` containing `funcmask_constrained.nii.gz`, is required unless `cfg.allowNoQualityMask=true` explicitly opts out. No resampling is done: image geometry must already agree. The signal validity rule currently requires all time points finite, nonzero, and temporally nonconstant. That is a conservative design choice, not a validated universal neuroimaging quality rule. Unsmoothed BOLD is preferred for anatomical/small ROIs, especially near DBS leads. The code reports raw voxel counts, quality/valid voxel counts, retained weight mass, weighted effective N, coverage fractions, and NaN time series for empty/low-coverage ROIs. User-selected thresholds have no validated scientific default (`minCoverageFraction=0`, `minValidVoxels=1`).

Output: `roi_mean_timeseries.csv` (time + ROI columns), `roi_coverage.csv` (maps columns to original files and integer labels/stack indices), and `roi_provenance.json` (input hashes, semantics, QC and settings). Writes require `cfg.writeOutputs=true` and a new destination with an existing parent. Source paths in JSON are currently absolute and not yet relocation-robust. Names are `roi_001`, etc.; meaningful atlas names can be joined using `roi_coverage.csv`. This is **descriptive ROI BOLD extraction**, not statistical inference and not an ROI contrast-beta table.

### Example: scan-specific unsmoothed BOLD and atlas

```matlab
addpath('/Users/keithdodd/GitHub/CICADA/roi','-begin');
cfg=struct();
cfg.boldFile='/ABS/PATH/TO/UNSMOOTHED_BOLD.nii.gz';
cfg.roiFiles={'/ABS/PATH/TO/LABEL_ATLAS.nii.gz', ...
              '/ABS/PATH/TO/WEIGHTED_STACK.nii.gz'};
cfg.roiTypes={'labels','weighted'};
cfg.qualityMaskFile='/ABS/PATH/TO/funcmask_constrained.nii.gz';
cfg.dryRun=true;
P=cicada_roi_timeseries(cfg);
disp(P.coverage);
cfg.dryRun=false;
cfg.writeOutputs=true;
cfg.outputDir='/ABS/PATH/TO/NEW_RESULT_DIRECTORY';
R=cicada_roi_timeseries(cfg);
```

Source core: `roi/helper_functions/cicada_roi_timeseries_v02.m` is byte-identical to the user-tested v0.2 candidate. The public entry point is new and requires the bounded smoke test before the integration is considered locally verified. Synthetic v0.2 tests passed 8/8 on work Mac, and the user reports the laptop real-data binary/0-1-weighted test also passed. Neither establishes validation for real fractional probabilistic atlases, geometric distortion near leads, portability, or ROI statistical inference.

Do not add a Brainnetome dependency; arbitrary user-supplied atlases are the intended workflow. Keep ROI inputs and exported results outside CICADA-managed source and derivatives trees.
