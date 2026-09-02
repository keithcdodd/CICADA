#!/usr/bin/env python3
"""Public-safe static and numeric checks for output-smoothing geometry."""

from math import isfinite, log, sqrt
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[2]
HELPERS = REPO_ROOT / "helper_functions"
FWHM_PER_SIGMA = 2 * sqrt(2 * log(2))


def sigma_vox(fwhm_mm, voxel_sizes_mm):
    assert isfinite(fwhm_mm) and fwhm_mm >= 0
    assert len(voxel_sizes_mm) == 3
    assert all(isfinite(v) and v > 0 for v in voxel_sizes_mm)
    sigma_mm = fwhm_mm / FWHM_PER_SIGMA
    return [sigma_mm / v for v in voxel_sizes_mm]


def smoothing_tag(fwhm_mm):
    assert isfinite(fwhm_mm) and fwhm_mm >= 0
    value_text = format(fwhm_mm, ".15g")
    return "s" + value_text.replace(".", "p").replace("+", "").replace("-", "m")


def main():
    helper = (HELPERS / "cicada_fwhm_mm_to_sigma_vox.m").read_text()
    tag_helper = (HELPERS / "cicada_smoothing_fwhm_tag.m").read_text()
    detrend = (HELPERS / "detrend_filter_smooth.m").read_text()
    network = (HELPERS / "network_identifiability.m").read_text()
    seed = (HELPERS / "seed_based_conn.m").read_text()
    group = (REPO_ROOT / "group_qc" / "cicada_group_qc.m").read_text()
    end_to_end = (
        REPO_ROOT / "validation" / "scripts" /
        "test_output_smoothing_end_to_end.m"
    ).read_text()

    assert "2 * sqrt(2 * log(2))" in helper
    assert "sigma_mm ./ voxel_sizes_mm" in helper
    assert "'numel',3,'finite','positive'" in helper
    assert "sprintf('%.15g'" in tag_helper
    assert "strrep(value_text, '.', 'p')" in tag_helper
    for active_path in (detrend, network, seed):
        assert "cicada_fwhm_mm_to_sigma_vox" in active_path
        assert "imgaussfilt3" in active_path
        assert "sigma_vox" in active_path
    assert "round(mean(first_image_nifi_info.PixelDimensions(1:3)))" not in group
    assert "isscalar(smoothing_kernel)" in group
    assert "~isnan(smoothing_kernel)" in group
    assert "~isinf(smoothing_kernel)" in group
    assert "num2str(round(smoothing_kernel))" not in detrend
    assert "cicada_smoothing_fwhm_tag" in detrend
    assert "detrend_filter_smooth(" in end_to_end
    assert "niftiwrite(" in end_to_end
    assert "niftiread(" in end_to_end
    assert "imgaussfilt3(" in end_to_end
    assert smoothing_tag(6.0) == "s6"
    assert smoothing_tag(6.4) == "s6p4"
    assert smoothing_tag(6.4) != smoothing_tag(6.49)

    # 2-mm legacy equivalence: exact FWHM conversion differs from the old
    # rounded 2.355 constant by less than 0.0001 voxel sigma.
    isotropic = sigma_vox(6.0, [2.0, 2.0, 2.0])
    legacy = (6.0 / 2.0) / 2.355
    assert all(abs(value - legacy) < 1e-4 for value in isotropic)

    # Anisotropic and non-integer voxels recover the same physical FWHM on
    # every axis and retain unequal voxel-space sigmas where required.
    for fwhm, voxels in ((6.0, [2.0, 2.0, 3.0]),
                         (5.5, [2.4, 2.4, 2.7])):
        sigmas = sigma_vox(fwhm, voxels)
        recovered = [s * v * FWHM_PER_SIGMA
                     for s, v in zip(sigmas, voxels)]
        assert all(abs(value - fwhm) < 1e-12 for value in recovered)
        assert sigmas[0] != sigmas[2]

    for invalid_fwhm, invalid_voxels in (
        (-1.0, [2.0, 2.0, 2.0]),
        (float("nan"), [2.0, 2.0, 2.0]),
        (6.0, [2.0, 0.0, 2.0]),
        (6.0, [2.0, float("inf"), 2.0]),
        (6.0, [2.0, 2.0]),
    ):
        try:
            sigma_vox(invalid_fwhm, invalid_voxels)
        except AssertionError:
            pass
        else:
            raise AssertionError("Invalid smoothing geometry was accepted")

    print("PASS: output-smoothing geometry configuration and numeric checks.")


if __name__ == "__main__":
    main()
