# Physical output-smoothing geometry validation

CICADA accepts requested Gaussian smoothing widths as physical FWHM in
millimeters. The active MATLAB output/QC smoothing paths now convert that
single physical width to a three-element voxel-space sigma:

```text
sigma_mm = FWHM_mm / (2 * sqrt(2 * log(2)))
sigma_vox = sigma_mm ./ [voxel_x_mm voxel_y_mm voxel_z_mm]
```

This preserves the requested physical width independently on x, y, and z for
isotropic, anisotropic, and non-integer voxel dimensions. Inputs are rejected
unless FWHM is finite and nonnegative and all three spatial voxel dimensions
are finite and positive. Runtime output records requested FWHM, voxel sizes,
and resulting voxel-space sigmas.

The validation suite covers:

- numerical equivalence to the historical scalar calculation for ordinary
  2 x 2 x 2 mm data (within 0.0001 voxel sigma, accounting for replacement of
  the rounded constant 2.355 by the exact conversion);
- recovery of the same physical FWHM on every axis for 2 x 2 x 3 mm data;
- preservation of non-integer 2.4 x 2.4 x 2.7 mm voxel sizes;
- rejection of negative/non-finite FWHM and invalid voxel dimensions; and
- backward-compatible integer output tags plus distinct non-rounded tags for
  non-integer requested kernels;
- end-to-end synthetic NIfTI execution through `detrend_filter_smooth`,
  including expected output values, voxel metadata, and filenames; and
- an explicit zero-smoothing case that confirms no Gaussian filter is applied.

Run the public-safe configuration/numeric check from the repository root:

```bash
python validation/scripts/test_output_smoothing_configuration.py
```

When MATLAB with Image Processing Toolbox is available, also run the native
helper and end-to-end tests:

```matlab
test_output_smoothing_geometry
test_output_smoothing_end_to_end
```

The end-to-end test creates only temporary synthetic impulse images and removes
them automatically. It does not access participant data.

The explicitly retained `network_identifiability_old.m` legacy implementation
is not changed by this hardening update.
