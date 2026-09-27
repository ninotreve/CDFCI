#!/usr/bin/env python3
"""Summarize the large CDFCI trajectory and timing-pair benchmark."""

from __future__ import annotations

import argparse
import csv
import json
import math
import statistics
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

REFERENCE_RESOLUTION_FLOOR = 1e-7

def first_crossing(rows: list[dict], key: str, tolerance: float) -> dict | None:
    return next((row for row in rows if float(row[key]) <= tolerance), None)


def reference_resolved_xy(xs: list[float], errors: list[float]) -> tuple[list[float], list[float]]:
    pairs = [(x, error) for x, error in zip(xs, errors)
             if error >= REFERENCE_RESOLUTION_FLOOR]
    return [x for x, _ in pairs], [error for _, error in pairs]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("input", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument(
        "--reference-energy", type=float,
        help="Override the energy stored in the trajectory when reanalysing an existing run.",
    )
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    data = json.loads(args.input.read_text(encoding="utf-8"))
    stored_reference = float(data["reference_energy"])
    reference = (stored_reference if args.reference_energy is None
                 else args.reference_energy)
    all_rows = []
    for original in data["trajectory"]:
        row = dict(original)
        row["raw_error"] = abs(float(row["variational_energy"]) - reference)
        row["external_only_error"] = abs(
            float(row["variational_energy"])
            + float(row["external_correction"])
            - reference
        )
        row["corrected_error"] = abs(float(row["corrected_energy"]) - reference)
        all_rows.append(row)
    rows = [row for row in all_rows if row["status"] == "ok"]
    if not rows:
        raise RuntimeError("large-scale trajectory has no valid correction records")
    reference_run = data.get("reference_run") or {}
    reference_tail_change = reference_run.get("tail_energy_change")

    threshold_rows = []
    for tolerance in (1e-3, 1e-4, 1e-5, 1e-6):
        raw = first_crossing(rows, "raw_error", tolerance)
        corrected = first_crossing(rows, "corrected_error", tolerance)
        threshold_rows.append({
            "tolerance_hartree": tolerance,
            "raw_iteration": raw["iteration"] if raw else "",
            "corrected_iteration": corrected["iteration"] if corrected else "",
            "raw_hamiltonian_columns": raw["hamiltonian_columns"] if raw else "",
            "corrected_hamiltonian_columns": corrected["hamiltonian_columns"] if corrected else "",
            "raw_wall_seconds": raw["raw_wall_seconds"] if raw else "",
            "corrected_wall_seconds": corrected["corrected_wall_seconds"] if corrected else "",
            "iteration_speedup": (float(raw["iteration"]) / float(corrected["iteration"]))
                if raw and corrected else "",
            "wall_time_speedup": (float(raw["raw_wall_seconds"]) /
                                  float(corrected["corrected_wall_seconds"]))
                if raw and corrected and float(corrected["corrected_wall_seconds"]) > 0 else "",
            "reference_tail_below_tenth_tolerance": (
                float(reference_tail_change) <= tolerance / 10
                if reference_tail_change is not None else "unknown"),
        })
    table_path = args.output / "05_target_tolerances.csv"
    with table_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(threshold_rows[0]))
        writer.writeheader()
        writer.writerows(threshold_rows)

    trajectory_fields = (
        "iteration", "stored_determinants", "stored_residual_entries",
        "stored_wavefunction_entries",
        "variational_energy", "corrected_energy", "raw_error", "corrected_error",
        "raw_wall_seconds", "corrected_wall_seconds", "correction_seconds",
        "interval_total_seconds", "interval_correction_overhead",
        "cumulative_correction_overhead", "hamiltonian_columns", "status",
    )
    with (args.output / "05_trajectory.csv").open(
            "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=trajectory_fields,
                                extrasaction="ignore")
        writer.writeheader()
        writer.writerows(all_rows)

    overheads = []
    for pair in data["timings"]:
        baseline = float(pair["baseline"]["seconds"])
        enabled = float(pair["enabled"]["seconds"])
        overheads.append((enabled - baseline) / baseline)
    run = data["trajectory_run"]
    correction_seconds = float(run["correction_seconds"])
    parent_seconds = float(run["seconds"]) - correction_seconds
    n_spatial = int(data["spin_orbitals"]) // 2
    n_alpha = (int(data["electrons"]) + int(data.get("ms2", 0))) // 2
    n_beta = int(data["electrons"]) - n_alpha
    determinant_dimension = math.comb(n_spatial, n_alpha) * math.comb(n_spatial, n_beta)
    paired_energy_differences = [
        abs(float(row["variational_energy"]) -
            float(row["baseline_variational_energy"])) for row in rows
    ]
    improved = [float(row["corrected_error"]) < float(row["raw_error"]) for row in rows]
    overshoots = [row for row in rows
                  if float(row["corrected_energy"]) < reference]
    summary = {
        "fcidump": data["fcidump"],
        "spatial_orbitals": n_spatial,
        "spin_orbitals": data["spin_orbitals"],
        "electrons": data["electrons"],
        "fixed_ms_determinant_dimension": str(determinant_dimension),
        "reference_energy": reference,
        "reference_kind": (data.get("reference_kind", "supplied")
                           if args.reference_energy is None else "analysis_override"),
        "trajectory_stored_reference_energy": stored_reference,
        "reference_iterations": reference_run.get("iterations"),
        "reference_tail_energy_change": reference_tail_change,
        "trajectory_points": len(rows),
        "plot_reference_resolution_floor_hartree": REFERENCE_RESOLUTION_FLOOR,
        "invalid_correction_points": len(all_rows) - len(rows),
        "improved_point_fraction": sum(improved) / len(improved),
        "overshoot_points": len(overshoots),
        "first_overshoot_iteration": overshoots[0]["iteration"] if overshoots else None,
        "median_relative_runtime_overhead": statistics.median(overheads) if overheads else None,
        "paired_relative_runtime_overheads": overheads,
        "direct_correction_seconds": correction_seconds,
        "direct_relative_runtime_overhead": (
            correction_seconds / parent_seconds if parent_seconds > 0 else None),
        "max_paired_variational_energy_difference": max(paired_energy_differences),
        "target_tolerances": threshold_rows,
    }
    (args.output / "05_large_scale_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8")

    fig, axes = plt.subplots(1, 3, figsize=(15.6, 4.4))
    raw_error = [float(row["raw_error"]) for row in rows]
    corrected_error = [float(row["corrected_error"]) for row in rows]
    stored = [int(row["stored_determinants"]) for row in rows]
    wavefunction_size = [int(row["stored_wavefunction_entries"]) for row in rows]
    raw_stored, raw_stored_error = reference_resolved_xy(stored, raw_error)
    corrected_stored, corrected_stored_error = reference_resolved_xy(stored, corrected_error)
    raw_time, raw_time_error = reference_resolved_xy(
        [float(row["raw_wall_seconds"]) for row in rows], raw_error)
    corrected_time, corrected_time_error = reference_resolved_xy(
        [float(row["corrected_wall_seconds"]) for row in rows], corrected_error)

    axes[0].loglog(raw_stored, raw_stored_error, "-", linewidth=1.5, label="raw CDFCI")
    axes[0].loglog(corrected_stored, corrected_stored_error, "-", linewidth=1.5,
                   label="corrected CDFCI")
    axes[0].set_title("A. Error vs stored determinants", loc="left")
    axes[0].set(xlabel="stored determinants", ylabel="absolute energy error (Ha)")
    axes[0].set_ylim(bottom=REFERENCE_RESOLUTION_FLOOR)
    axes[0].legend()

    axes[1].loglog(raw_time, raw_time_error, "-", linewidth=1.5, label="raw CDFCI")
    axes[1].loglog(corrected_time, corrected_time_error, "-", linewidth=1.5,
                   label="corrected CDFCI")
    axes[1].set_title("B. Error vs wall time", loc="left")
    axes[1].set(xlabel="wall time (s)", ylabel="absolute energy error (Ha)")
    axes[1].set_ylim(bottom=REFERENCE_RESOLUTION_FLOOR)
    axes[1].legend()

    axes[2].semilogx(
        wavefunction_size,
        [float(row["interval_correction_overhead"]) for row in rows],
        "-", linewidth=1.5, color="tab:green")
    axes[2].set_title("C. Correction overhead", loc="left")
    axes[2].set(
        xlabel="stored wavefunction entries",
        ylabel="correction time / interval total time",
    )
    for ax in axes:
        ax.grid(True, which="both", alpha=0.25)
    fig.suptitle(
        f"Reference energy = {reference:.10f} Ha; "
        r"only $|\Delta E| \geq 10^{-7}$ Ha shown",
        fontsize=11,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    fig.savefig(args.output / "05_correction_advantage.pdf")
    fig.savefig(args.output / "05_correction_advantage.png", dpi=220)
    fig.savefig(args.output / "05_correction_advantage_reference_resolved.pdf")
    fig.savefig(args.output / "05_correction_advantage_reference_resolved.png", dpi=220)
    plt.close(fig)
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
