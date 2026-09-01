#!/usr/bin/env python3

import argparse
import csv
import sys
import statistics
from pathlib import Path


def parse_args():
    p = argparse.ArgumentParser(
        description=(
            "Reconstruct CICADA classifier-validation metrics from saved "
            "IC_auto_checker.csv and IC_manual_checker.csv files."
        )
    )
    p.add_argument(
        "--master",
        required=True,
        type=Path,
        help="Path to cicada_runs_master.csv",
    )
    p.add_argument(
        "--data-root",
        required=True,
        type=Path,
        help="Root CICADA derivatives directory containing sub-* folders",
    )
    p.add_argument(
        "--output-dir",
        required=True,
        type=Path,
    )
    p.add_argument(
        "--roles",
        nargs="+",
        default=["manual"],
        help=("Roles to include from master CSV. Historical roles supported: "
              "manual (consensus labels) and fixtraining (KD individual labels). "
              "Default: manual"),
    )
    p.add_argument(
        "--allow-output-inside-git",
        action="store_true",
        help=(
            "DANGEROUS OVERRIDE: allow participant-derived validation outputs "
            "to be written inside a Git repository. Disabled by default."
        ),
    )
    return p.parse_args()


def normalize_ic(value):
    """Normalize PotentialICs to an integer IC number."""
    return int(float(str(value).strip()))


def normalize_label(value, filename, ic):
    try:
        label = int(float(str(value).strip()))
    except Exception as exc:
        raise ValueError(
            f"Invalid SignalLabel '{value}' for IC {ic} in {filename}"
        ) from exc

    if label not in (0, 1):
        raise ValueError(
            f"SignalLabel must be 0 or 1; got {label} "
            f"for IC {ic} in {filename}"
        )

    return label


def read_checker(path):
    """
    Read an IC checker CSV and return a dictionary keyed by PotentialICs.

    Automatic and manual checker rows are not assumed to have the same
    ordering. Pairing is always performed by IC number.
    """
    rows = {}

    with path.open(newline="", encoding="utf-8-sig") as f:
        reader = csv.DictReader(f)

        required = {"PotentialICs", "SignalLabel"}
        missing = required - set(reader.fieldnames or [])

        if missing:
            raise ValueError(
                f"{path} is missing required columns: {sorted(missing)}"
            )

        for row in reader:
            ic = normalize_ic(row["PotentialICs"])

            if ic in rows:
                raise ValueError(
                    f"Duplicate PotentialICs={ic} in {path}"
                )

            rows[ic] = {
                "SignalLabel": normalize_label(
                    row["SignalLabel"], path, ic
                ),
                "raw": row,
            }

    if not rows:
        raise ValueError(f"No IC rows found in {path}")

    return rows


def safe_div(num, den):
    return num / den if den else None


def mean_sd(values):
    """Return mean and sample SD after dropping None values."""
    clean = [float(v) for v in values if v is not None]
    if not clean:
        return None, None
    mean_value = statistics.mean(clean)
    sd_value = statistics.stdev(clean) if len(clean) > 1 else 0.0
    return mean_value, sd_value


def calculate_metrics(rows):
    """
    Manual SignalLabel is treated as the gold standard.

    Positive class = signal.
    Negative class = noise.
    """
    tp = sum(
        1 for r in rows
        if r["manual_label"] == 1 and r["auto_label"] == 1
    )
    fn = sum(
        1 for r in rows
        if r["manual_label"] == 1 and r["auto_label"] == 0
    )
    fp = sum(
        1 for r in rows
        if r["manual_label"] == 0 and r["auto_label"] == 1
    )
    tn = sum(
        1 for r in rows
        if r["manual_label"] == 0 and r["auto_label"] == 0
    )

    signal_sensitivity = safe_div(tp, tp + fn)
    signal_ppv = safe_div(tp, tp + fp)

    if signal_sensitivity is None or signal_ppv is None:
        signal_f1 = None
    elif signal_sensitivity + signal_ppv == 0:
        signal_f1 = 0.0
    else:
        signal_f1 = (
            2 * signal_sensitivity * signal_ppv
            / (signal_sensitivity + signal_ppv)
        )

    return {
        "n_ics": tp + fn + fp + tn,
        "manual_signal": tp + fn,
        "manual_noise": tn + fp,
        "TP": tp,
        "FN": fn,
        "FP": fp,
        "TN": tn,
        "signal_sensitivity": signal_sensitivity,
        "signal_ppv": signal_ppv,
        "signal_f1": signal_f1,
        "noise_sensitivity": safe_div(tn, tn + fp),
        "npv": safe_div(tn, tn + fn),
        "accuracy": safe_div(tp + tn, tp + tn + fp + fn),
    }


def write_csv(path, rows, fieldnames):
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)



def find_git_root(path):
    """
    Return the nearest ancestor that appears to be a Git work tree.

    A .git entry may be either a directory (ordinary clone) or a file
    (worktree/submodule), so existence is sufficient for this safety check.
    """
    path = path.resolve()

    for candidate in (path, *path.parents):
        if (candidate / ".git").exists():
            return candidate

    return None


def gold_standard_file(scan_dir, role):
    """
    Return the historically appropriate manual-label file and provenance.

    manual:
        Final consensus classification used for the held-out validation
        subset. Requires IC_manual_checker.csv.

    fixtraining:
        Historical FIX-training labels were KD's individual manual labels.
        Naming was not fully standardized at the time:
          1) Prefer IC_manual_kd_checker.csv when present.
          2) Otherwise use IC_manual_checker.csv as the legacy/default
             filename for KD's individual labels.

        The fallback is explicitly recorded in output provenance so it is
        never confused with the later held-out consensus convention.
    """
    if role == "manual":
        path = scan_dir / "IC_manual_checker.csv"
        return path, "consensus_KD_LS_MM"

    if role == "fixtraining":
        kd_path = scan_dir / "IC_manual_kd_checker.csv"
        legacy_path = scan_dir / "IC_manual_checker.csv"

        if kd_path.is_file():
            return kd_path, "KD_individual_explicit_filename"

        if legacy_path.is_file():
            return (
                legacy_path,
                "KD_individual_legacy_default_filename",
            )

        # Return the preferred explicit filename so the missing-file report
        # points to the intended historical label type.
        return kd_path, "KD_individual_missing"

    raise ValueError(
        f"Unsupported validation role '{role}'. "
        "Use role 'manual' or 'fixtraining' for historical reconstruction."
    )


def main():
    args = parse_args()

    args.master = args.master.expanduser().resolve()
    args.data_root = args.data_root.expanduser().resolve()
    args.output_dir = args.output_dir.expanduser().resolve()

    wanted_roles = set(args.roles)

    # ------------------------------------------------------------------
    # Privacy / data-governance guard
    #
    # Historical "manual" and "fixtraining" runs contain participant-level
    # and/or IC-level information derived from human-subject research data.
    # By default, refuse to write those outputs anywhere inside a Git work
    # tree. This prevents an accidental `git add` / public push.
    # ------------------------------------------------------------------
    human_subject_roles = {"manual", "fixtraining"}
    git_root = find_git_root(args.output_dir)

    if (
        wanted_roles & human_subject_roles
        and git_root is not None
        and not args.allow_output_inside_git
    ):
        raise RuntimeError(
            "\nPRIVACY SAFETY STOP\n"
            "Participant-derived CICADA validation output was requested "
            "inside a Git repository.\n\n"
            f"Requested output: {args.output_dir}\n"
            f"Git repository:    {git_root}\n\n"
            "Write these results to a private directory outside Git, for "
            "example:\n"
            "  $HOME/CICADA_private_validation/<date>/<benchmark>\n\n"
            "The override --allow-output-inside-git exists only for an "
            "intentional, reviewed use case."
        )

    args.output_dir.mkdir(parents=True, exist_ok=True)

    # utf-8-sig safely handles a UTF-8 BOM on the first CSV header.
    with args.master.open(newline="", encoding="utf-8-sig") as f:
        master_rows = list(csv.DictReader(f))

    if not master_rows:
        raise RuntimeError(
            f"No rows found in master CSV: {args.master}"
        )

    required_master_cols = {
        "sub_id",
        "ses_id",
        "task_name",
        "role",
        "exclude",
    }

    actual_master_cols = set(master_rows[0].keys())
    missing_master_cols = required_master_cols - actual_master_cols

    if missing_master_cols:
        raise RuntimeError(
            "Master CSV is missing required columns: "
            f"{sorted(missing_master_cols)}\n"
            f"Columns actually found: {list(master_rows[0].keys())}"
        )

    selected = [
        r for r in master_rows
        if r.get("role", "").strip() in wanted_roles
        and r.get("exclude", "0").strip() != "1"
        and r.get("task_name", "").strip()
        in {"rest", "foodpics_run-01"}
    ]

    print(f"Selected {len(selected)} scans")
    print(f"Roles: {sorted(wanted_roles)}")

    if wanted_roles == {"manual"}:
        print("Gold labels: final consensus IC_manual_checker.csv")
    elif wanted_roles == {"fixtraining"}:
        print(
            "Gold labels: KD individual labels; prefer "
            "IC_manual_kd_checker.csv, otherwise use legacy "
            "IC_manual_checker.csv naming"
        )
    else:
        print("Gold labels: role-specific historical manual files")

    ic_rows = []
    missing_files = []

    for m in selected:
        sub = str(m["sub_id"]).strip()
        ses = int(float(m["ses_id"]))
        task = m["task_name"].strip()
        role = m["role"].strip()

        scan_dir = (
            args.data_root
            / f"sub-{sub}"
            / f"ses-{ses:02d}"
            / task
            / "ic_auto_selection"
        )

        auto_file = scan_dir / "IC_auto_checker.csv"
        manual_file, gold_label_source = gold_standard_file(
            scan_dir, role
        )

        missing = [
            str(p)
            for p in (auto_file, manual_file)
            if not p.is_file()
        ]

        if missing:
            missing_files.append(
                {
                    "sub_id": sub,
                    "ses_id": ses,
                    "task_name": task,
                    "role": role,
                    "gold_label_source": gold_label_source,
                    "missing": "; ".join(missing),
                }
            )
            continue

        auto = read_checker(auto_file)
        manual = read_checker(manual_file)

        auto_ics = set(auto)
        manual_ics = set(manual)

        if auto_ics != manual_ics:
            auto_only = sorted(auto_ics - manual_ics)
            manual_only = sorted(manual_ics - auto_ics)

            raise RuntimeError(
                "\nIC-set mismatch:\n"
                f"  subject: sub-{sub}\n"
                f"  task: {task}\n"
                f"  auto only: {auto_only}\n"
                f"  manual only: {manual_only}\n"
            )

        # Sort only for stable output. Pairing itself is by IC number.
        for ic in sorted(auto_ics):
            auto_label = auto[ic]["SignalLabel"]
            manual_label = manual[ic]["SignalLabel"]

            if manual_label == 1 and auto_label == 1:
                outcome = "TP"
            elif manual_label == 1 and auto_label == 0:
                outcome = "FN"
            elif manual_label == 0 and auto_label == 1:
                outcome = "FP"
            else:
                outcome = "TN"

            ic_rows.append(
                {
                    "sub_id": sub,
                    "ses_id": ses,
                    "task_name": task,
                    "role": role,
                    "gold_label_source": gold_label_source,
                    "PotentialICs": ic,
                    "auto_label": auto_label,
                    "manual_label": manual_label,
                    "discordant": int(auto_label != manual_label),
                    "outcome": outcome,
                    "auto_file": str(auto_file),
                    "manual_file": str(manual_file),
                }
            )

    if missing_files:
        missing_out = (
            args.output_dir / "missing_validation_files.csv"
        )
        write_csv(
            missing_out,
            missing_files,
            [
                "sub_id",
                "ses_id",
                "task_name",
                "role",
                "gold_label_source",
                "missing",
            ],
        )

        print("\nERROR: required validation files are missing.")
        print(f"See: {missing_out}")
        print("No historical metrics were written.")
        return 2

    rest_scans = {
        (r["sub_id"], r["ses_id"])
        for r in ic_rows
        if r["task_name"] == "rest"
    }

    task_scans = {
        (r["sub_id"], r["ses_id"])
        for r in ic_rows
        if r["task_name"] == "foodpics_run-01"
    }

    print(f"Usable rest scans: {len(rest_scans)}")
    print(f"Usable task scans: {len(task_scans)}")

    source_counts = {}
    scan_sources_seen = set()

    for r in ic_rows:
        scan_key = (
            r["sub_id"],
            r["ses_id"],
            r["task_name"],
            r["gold_label_source"],
        )
        scan_sources_seen.add(scan_key)

    for _, _, _, source in scan_sources_seen:
        source_counts[source] = source_counts.get(source, 0) + 1

    print("Gold-label provenance by scan:")
    for source in sorted(source_counts):
        print(f"  {source}: {source_counts[source]}")

    ic_out = (
        args.output_dir / "historical_ic_level_results.csv"
    )

    write_csv(
        ic_out,
        ic_rows,
        [
            "sub_id",
            "ses_id",
            "task_name",
            "role",
            "gold_label_source",
            "PotentialICs",
            "auto_label",
            "manual_label",
            "discordant",
            "outcome",
            "auto_file",
            "manual_file",
        ],
    )

    scan_keys = sorted(
        {
            (
                r["sub_id"],
                r["ses_id"],
                r["task_name"],
                r["role"],
                r["gold_label_source"],
            )
            for r in ic_rows
        }
    )

    scan_metrics = []

    for sub, ses, task, role, gold_label_source in scan_keys:
        rows = [
            r for r in ic_rows
            if (
                r["sub_id"] == sub
                and r["ses_id"] == ses
                and r["task_name"] == task
                and r["role"] == role
                and r["gold_label_source"] == gold_label_source
            )
        ]

        metrics = calculate_metrics(rows)

        scan_metrics.append(
            {
                "sub_id": sub,
                "ses_id": ses,
                "task_name": task,
                "role": role,
                "gold_label_source": gold_label_source,
                **metrics,
            }
        )

    metric_fields = [
        "n_ics",
        "manual_signal",
        "manual_noise",
        "TP",
        "FN",
        "FP",
        "TN",
        "signal_sensitivity",
        "signal_ppv",
        "signal_f1",
        "noise_sensitivity",
        "npv",
        "accuracy",
    ]

    scan_out = (
        args.output_dir / "historical_scan_metrics.csv"
    )

    write_csv(
        scan_out,
        scan_metrics,
        [
            "sub_id",
            "ses_id",
            "task_name",
            "role",
            "gold_label_source",
        ]
        + metric_fields,
    )

    # ------------------------------------------------------------------
    # Dataset summaries
    #
    # "micro" metrics pool all ICs before constructing the confusion
    # matrix. These are useful regression summaries, but they are not the
    # same statistic reported in the CICADA publication.
    #
    # "mean_subject" metrics first compute each scan's confusion-matrix
    # metrics and then average those metrics across scans. The publication
    # reports mean subject accuracies, so these are the primary historical
    # reproduction statistics.
    # ------------------------------------------------------------------

    micro_summary_rows = []

    groups = [
        (
            "rest",
            [
                r for r in ic_rows
                if r["task_name"] == "rest"
            ],
        ),
        (
            "foodpics_run-01",
            [
                r for r in ic_rows
                if r["task_name"] == "foodpics_run-01"
            ],
        ),
        ("combined_descriptive", ic_rows),
    ]

    for name, rows in groups:
        metrics = calculate_metrics(rows)

        n_scans = len(
            {
                (
                    r["sub_id"],
                    r["ses_id"],
                    r["task_name"],
                )
                for r in rows
            }
        )

        micro_summary_rows.append(
            {
                "dataset": name,
                "summary_type": "pooled_IC_micro",
                "n_scans": n_scans,
                **metrics,
            }
        )

    micro_summary_out = (
        args.output_dir / "historical_summary_metrics_pooled_IC.csv"
    )

    write_csv(
        micro_summary_out,
        micro_summary_rows,
        ["dataset", "summary_type", "n_scans"] + metric_fields,
    )

    subject_metric_names = [
        "signal_sensitivity",
        "signal_ppv",
        "signal_f1",
        "noise_sensitivity",
        "npv",
        "accuracy",
    ]

    subject_summary_rows = []

    subject_groups = [
        (
            "rest",
            [
                r for r in scan_metrics
                if r["task_name"] == "rest"
            ],
        ),
        (
            "foodpics_run-01",
            [
                r for r in scan_metrics
                if r["task_name"] == "foodpics_run-01"
            ],
        ),
        # Descriptive only. Rest and task scans often come from the same
        # participants, so do not treat this 60-scan row as 60 independent
        # participants for inferential statistics.
        ("combined_descriptive", scan_metrics),
    ]

    for name, rows in subject_groups:
        out = {
            "dataset": name,
            "summary_type": "mean_subject_macro",
            "n_scans": len(rows),
        }

        for metric in subject_metric_names:
            mean_value, sd_value = mean_sd(
                [r[metric] for r in rows]
            )
            out[f"{metric}_mean"] = mean_value
            out[f"{metric}_sd"] = sd_value

        subject_summary_rows.append(out)

    subject_summary_fields = [
        "dataset",
        "summary_type",
        "n_scans",
    ]

    for metric in subject_metric_names:
        subject_summary_fields.extend(
            [f"{metric}_mean", f"{metric}_sd"]
        )

    subject_summary_out = (
        args.output_dir / "historical_summary_metrics_mean_subject.csv"
    )

    write_csv(
        subject_summary_out,
        subject_summary_rows,
        subject_summary_fields,
    )

    # Path-sanitized IC-level reference. Still PRIVATE human-subject-derived data.
    private_ic_out = (
        args.output_dir / "historical_ic_level_reference_private.csv"
    )

    write_csv(
        private_ic_out,
        [
            {
                k: r[k]
                for k in [
                    "sub_id",
                    "ses_id",
                    "task_name",
                    "role",
                    "gold_label_source",
                    "PotentialICs",
                    "auto_label",
                    "manual_label",
                    "discordant",
                    "outcome",
                ]
            }
            for r in ic_rows
        ],
        [
            "sub_id",
            "ses_id",
            "task_name",
            "role",
            "gold_label_source",
            "PotentialICs",
            "auto_label",
            "manual_label",
            "discordant",
            "outcome",
        ],
    )

    print("\nMean subject-level validation summary")
    print("=" * 100)

    for row in subject_summary_rows:
        print(
            f"{row['dataset']:22s} "
            f"scans={row['n_scans']:2d} "
            f"SS={row['signal_sensitivity_mean']:.4f}"
            f"±{row['signal_sensitivity_sd']:.4f} "
            f"SPV={row['signal_ppv_mean']:.4f}"
            f"±{row['signal_ppv_sd']:.4f} "
            f"F1={row['signal_f1_mean']:.4f}"
            f"±{row['signal_f1_sd']:.4f} "
            f"NS={row['noise_sensitivity_mean']:.4f}"
            f"±{row['noise_sensitivity_sd']:.4f} "
            f"NPV={row['npv_mean']:.4f}"
            f"±{row['npv_sd']:.4f} "
            f"OA={row['accuracy_mean']:.4f}"
            f"±{row['accuracy_sd']:.4f}"
        )

    print("\nPooled IC-level validation summary")
    print("=" * 100)

    for row in micro_summary_rows:
        print(
            f"{row['dataset']:22s} "
            f"scans={row['n_scans']:2d} "
            f"ICs={row['n_ics']:4d} "
            f"SS={row['signal_sensitivity']:.4f} "
            f"SPV={row['signal_ppv']:.4f} "
            f"F1={row['signal_f1']:.4f} "
            f"NS={row['noise_sensitivity']:.4f} "
            f"NPV={row['npv']:.4f} "
            f"OA={row['accuracy']:.4f}"
        )

    # For the historical held-out manual benchmark, verify the expected
    # accessible schizophrenia scan counts. This is a regression guard,
    # not a general requirement for other roles.
    if wanted_roles == {"manual"}:
        if len(rest_scans) != 30 or len(task_scans) != 30:
            raise RuntimeError(
                "Historical consensus manual benchmark expected 30 rest "
                "and 30 foodpics_run-01 scans, but found "
                f"{len(rest_scans)} rest and {len(task_scans)} task."
            )

    if wanted_roles == {"fixtraining"}:
        if len(rest_scans) != 10 or len(task_scans) != 10:
            raise RuntimeError(
                "Historical FIX-training benchmark expected 10 rest and "
                "10 foodpics_run-01 scans, but found "
                f"{len(rest_scans)} rest and {len(task_scans)} task."
            )

    discordant = sum(
        r["discordant"] for r in ic_rows
    )

    print(f"\nTotal discordant ICs: {discordant}")

    print("\nWrote:")
    print(f"  {ic_out}")
    print(f"  {private_ic_out}")
    print(f"  {scan_out}")
    print(f"  {subject_summary_out}")
    print(f"  {micro_summary_out}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
