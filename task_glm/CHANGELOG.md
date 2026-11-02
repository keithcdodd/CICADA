# Option 1 companion — change log

## v0.3.2-rc1 (unreleased, 2026-11-01)
- Single provenance version-string correction to `0.3.2-rc1-unreleased`.
- Updated README, validation summary, release checklist and GPLv3 distribution license.
- No REML, GLS, preflight, NIfTI export, or safety-algorithm changes.
- Bounded four-voxel read-only MATLAB version/entry-point check **passed**. Following that check, release-candidate documentation was revised to reflect identical Auto/Manual IC indices and the foodpics-only task inventory. Source `.m` files were **not changed** by that documentation review.

## v0.3.1-dev
- Fixed SHA-256 file reader spurious EOF error; fixture passed in user's MATLAB.
- 16-voxel fingerprinted preview and 11-file export passed.

## v0.3.0-dev
- Added SHA-256 input/code checksums to new fits and verified exports.
- Numerical fitting methods unchanged.

## v0.2-dev
- Separated public `task_glm.m` / export entry points from internal helper functions.
- 32-voxel three-contrast interface test passed.
