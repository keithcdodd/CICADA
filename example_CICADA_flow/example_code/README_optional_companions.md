# Optional CICADA companion examples (experimental)

Most CICADA workflows start with `example_Auto_CICADA.m` or `example_fmriprep_auto_CICADA.m`, optionally proceed through Manual CICADA, and use `example_Group_CICADA.m` for group-level processing/QC. Two **independent, optional** companion entry points are now available:

| What you want | Ready-to-edit example | MATLAB public function |
| --- | --- | --- |
| Task contrast from a newly verified, unsmoothed Group-prepared scan | `../../task_glm/examples/example_cicada_task_glm_group_matched.m` (relative to this folder) | `task_glm/cicada_task_glm.m` |
| Task contrast from the original scan with selected ICA noise modeled jointly | `../../task_glm/examples/example_cicada_task_glm_single.m` | `task_glm/cicada_task_glm.m` |
| Descriptive mean BOLD time series for a labeled atlas, binary masks, or weighted/4-D ROIs | `example_ROI_CICADA.m` | `roi/cicada_roi_timeseries.m` |

**No extra wrapper is required:** both public functions already accept one configuration struct and manage their internal helper paths. Use `cfg.dryRun=true` before attempting a full operation; the examples default to read-only preflight.

Group-matched task GLM **requires a new producer-attested matched Group recipe** linked to the exact scan/IC classification, with second-order polynomial drift, no spatial smoothing, and no voxelwise scaling. It does *not* accept unattested historical Group derivatives or `legacy_sequential` outputs. Setting `cfg.allowGroupExport=true` permits guarded first-level COPE and provisional VARCOPE export after explicit opt-in; it does not establish calibrated t/z/p inference. The task-GLM companion does not run ICA, Group processing, or a group-level statistical model.

ROI extraction is atlas-independent and reports a time-by-ROI CSV, voxel/weight-coverage CSV, and JSON provenance. It requires strict BOLD/ROI geometry alignment. Use `funcmask_constrained.nii.gz` as the canonical analysis-quality mask for new unified runs, or provide an explicitly aligned mask; legacy masks and MRI quality should still be assessed. For small anatomical ROIs prefer unsmoothed BOLD. ROI output is descriptive, **not** a contrast-beta or inferential result.

Add only `task_glm/` or `roi/` to MATLAB path for their respective public entry points; do not permanently add their `helper_functions/` folders. See `task_glm/README.md` and `roi/README.md` for detailed input contracts and limitations.
