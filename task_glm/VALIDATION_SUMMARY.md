# Option 1 companion — validation record as of 2026-11-01

**Scope:** original-time, joint task+drift+selected-ICA-noise design, external REML/FGLS companion. These are user-supplied local MATLAB results and earlier independently inspected numerical tests. They do not establish unconditional inferential validity.

## Verified computational behavior

- Exact reconstruction of FILM design: max absolute discrepancy **8.88e-16** (single-contrast sub-132 test).
- AR(1)/AR(2) cached vs reference engine on 64 matched voxels in each of sub-132 and sub-135: identical REML scores; max COPE difference **7.53e-13**; max relative provisional VARCOPE difference **4.25e-14**. Independent direct GLS, score and AR2 grid optimum checks passed.
- Sub-132 full-mask AR(1): **254379/254379** valid, **6.67 s** (user hardware); full NIfTI export **7** outputs, re-read verified.
- Sub-135 full-mask AR(1): **255526/255528** valid, two constant time-series voxels correctly excluded, **17.71 s**; full NIfTI export **7** outputs, re-read verified, attempted/valid masks checked.
- Auto and Manual preflight branches exercised. Manual sub-132 N=172, 52 noise ICs, X0=172x61, residual dimension=111; 64/64 bounded AR(1) fits passed. **The compared Auto and Manual noise-IC index sets were identical** (52 each; no differing-label run has been tested).
- v0.2 companion entry point and independent export entry point passed 32-voxel, 3-contrast in-memory test, including H-L=(H-N)-(L-N) identity error **8.88e-16**, no writes, path restoration.
- v0.3.1 provenance: SHA-256 known fixture PASSED; 16/16 read-only fit generated six input and nine code fingerprints; fingerprinted sparse export **11 files** re-read verified. Initial EOF bug in SHA-256 reader corrected. v0.3.2-rc1 **passed its four-voxel read-only smoke test** on sub-135, with 4/4 valid, `0.3.2-rc1-unreleased` recorded, `PREVIEW_ONLY` confirmed, and no files written. The statistical algorithms did not change.
- Exporter safety: wrong subject, changed design, managed CICADA destination, and existing destination were all rejected in source-tested MATLAB calls.

- Local inventory found **three** `hrf_conditions_raw.csv` scan directories (sub-132/sub-135/sub-140), all `foodpics_run-01`; no different task or different TR was established from this inventory. No general-task validation is claimed.

## Statistical calibration: not an inferential release

- Fixed-ICA Gaussian AR(1) nominal 5% null tests, 2000 replicates per rho: rejection **5.20%, 5.45%, 5.45%, 6.00%** for true r1 = 0, 0.3, 0.6, 0.8 respectively.
- Fixed-ICA nonzero task effect simulations generally recovered known coefficient with low bias; at true r1=0.8, nominal 95% CI coverage was **93.0% and 92.9%** at the two signal strengths (1000 replicates each).
- Paired oracle-vs-feasible at r1=0.4: **94.85% vs 93.85%** 95% CI coverage; at r1=0.8 **95.50% vs 93.95%**. Finite-sample uncertainty due to covariance estimation is unresolved.
- True AR(2) null simulations with r1=0.4, kappa2=+0.2: oracle/AR1/AR2 false rejection **4.20% / 7.85% / 5.90%**; kappa2=-0.2: **5.00% / 3.55% / 6.20%**.
- Stationary AR(1) vs changing temporal covariance stress test: with rho changing 0.1 -> 0.7, oracle/AR1/AR2 false rejection **4.75% / 7.85% / 7.10%**. Neither stationary covariance family is robustly calibrated under this test.
- Sub-132 real-data comparison with FILM over 254379 paired mask voxels: COPE Pearson **0.9691**, median REML/FILM VARCOPE ratio **1.294**; this is agreement, not validation of p-values.

**Release distinction:** exploratory COPE/provisional VARCOPE under explicitly bounded assumptions is viable; t/z/p/CI maps are deliberately excluded pending end-to-end adaptive ICA and covariance-calibration studies. These experiments have limited scenario coverage and do not amount to proof of generality.
