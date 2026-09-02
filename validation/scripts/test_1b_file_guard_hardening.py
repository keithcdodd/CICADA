#!/usr/bin/env python3
"""Public-safe static regression checks for CICADA roadmap item 1B."""

from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]


def main():
    qc3 = (ROOT / "basescripts" / "CICADA_3_QC.m").read_text()
    file_qc = (ROOT / "helper_functions" / "CICADA_fileQC.m").read_text()
    group_qc = (ROOT / "group_qc" / "cicada_group_qc.m").read_text()
    auto = (ROOT / "basescripts" / "CICADA_2_AutoLabeling.m").read_text()
    smoothness = (ROOT / "helper_functions" / "estimate_ic_smoothness.m").read_text()
    ident = (ROOT / "helper_functions" / "network_identifiability.m").read_text()

    # Direct CICADA_fileQC guards remain correct.
    assert "if ~isfile(denoised_file)" in file_qc
    assert "if ~isfile(compare_file)" in file_qc
    assert "if ~isfile(orig_file)" in file_qc

    # CICADA_3_QC validates cleaned input before dir() dereference.
    assert qc3.index("if ~isfile(cleaned_file)") < qc3.index(
        "cleaned_file_info = dir(cleaned_file)"
    )

    # Original-file lookup explicitly handles zero/multiple matches.
    assert "orig_file_info = dir(fullfile(cleaned_dir, '*orig*.nii.gz'))" in qc3
    assert "if isempty(orig_file_info)" in qc3
    assert "elseif numel(orig_file_info) > 1" in qc3
    assert qc3.index("if isempty(orig_file_info)") < qc3.index(
        "orig_file = fullfile(orig_file_info(1).folder"
    )

    # Default comparison is resolved by find_compare_file rather than an
    # unchecked wildcard dereference.
    assert "compare_file = '8p';" in qc3
    assert "compare_file_info = dir([cleaned_dir, '/*8p*'])" not in qc3
    assert "compare_file = find_compare_file(output_dir, compare_file, valid_tags);" in qc3

    # Group-QC comment reflects the actual poorly_improved expression.
    assert "FD, DVARS, & Spikiness" not in group_qc
    assert group_qc.count("% and decreased FD and DVARS") == 2

    # Auto-CICADA comment reflects the existing validated behavior.
    assert "not currently used as a marker" not in auto
    assert "absolute Spikiness value also exceeds 5" in auto
    assert "(Spikiness > 5)" in auto
    assert "'High_Spikiness'" in auto

    # Temporary-file names are unique and automatically cleaned.
    assert "ic_file_prefix = tempname;" in smoothness
    assert "onCleanup" in smoothness
    assert "ic_file_prefix = 'temp_ic';" not in smoothness

    assert "resampled_template = [tempname, '.nii.gz'];" in ident
    assert "onCleanup" in ident
    assert "fullfile(template_dir, 'resampled_template.nii.gz')" not in ident

    print(
        "PASS: CICADA 1B file-guard, documentation, "
        "and temporary-file hardening checks."
    )


if __name__ == "__main__":
    main()
