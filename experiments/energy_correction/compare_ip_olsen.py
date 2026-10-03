#!/usr/bin/env python3
"""Compare an IP-only baseline2 with an existing IP+Olsen run."""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


def compare(ip: dict, olsen: dict) -> tuple[dict, list[dict]]:
    if ip.get("benchmark_mode") != "ip-only":
        raise ValueError("baseline2 must be an ip-only benchmark")
    for data, expected in ((ip, False), (olsen, True)):
        correction = data["options"]["energy_correction"]
        if not correction["enabled"] or correction["olsen_enabled"] != expected:
            raise ValueError("comparison requires IP-only versus IP+Olsen")
    for key in ("fcidump", "electrons", "spin_orbitals", "ms2", "reference_energy"):
        if ip[key] != olsen[key]:
            raise ValueError(f"incompatible {key}")
    for key in ("num_iterations", "report_interval", "num_coordinates", "z_threshold",
                "max_wavefunction_size", "stopping_dx_threshold", "verbose"):
        if ip["options"][key] != olsen["options"][key]:
            raise ValueError(f"incompatible solver option {key}")
    ip_rows, olsen_rows = ip["trajectory"], olsen["trajectory"]
    if not ip_rows or len(ip_rows) != len(olsen_rows):
        raise ValueError("trajectories must have the same nonzero report count")
    rows = []
    for a, b in zip(ip_rows, olsen_rows):
        if a["iteration"] != b["iteration"] or a["hamiltonian_columns"] != b["hamiltonian_columns"]:
            raise ValueError("iteration/Hamiltonian-work mismatch")
        if a["internal_correction"] != 0:
            raise ValueError("baseline2 contains an Olsen correction")
        if a["status"] not in {"ok", "unchecked"} or b["status"] not in {"ok", "unchecked"}:
            raise ValueError("invalid correction status")
        ip_time, olsen_time = float(a["wall_seconds"]), float(b["corrected_wall_seconds"])
        if not all(math.isfinite(float(value)) for value in
                   (ip_time, olsen_time, a["variational_energy"], b["variational_energy"],
                    a["corrected_energy"], b["corrected_energy"])) or ip_time <= 0:
            raise ValueError("invalid time/energy record")
        rows.append({
            "iteration": a["iteration"],
            "ip_stored_determinants": a["stored_determinants"],
            "olsen_stored_determinants": b["stored_determinants"],
            "ip_stored_wavefunction_entries": a["stored_wavefunction_entries"],
            "olsen_stored_wavefunction_entries": b["stored_wavefunction_entries"],
            "ip_wall_seconds": ip_time,
            "ip_olsen_wall_seconds": olsen_time,
            "olsen_extra_wall_seconds": olsen_time - ip_time,
            "olsen_relative_overhead": (olsen_time - ip_time) / ip_time,
            "ip_report_seconds": a["correction_seconds"],
            "ip_olsen_report_seconds": b["correction_seconds"],
            "ip_variational_energy": a["variational_energy"],
            "ip_olsen_variational_energy": b["variational_energy"],
            "ip_corrected_energy": a["corrected_energy"],
            "ip_olsen_corrected_energy": b["corrected_energy"],
        })
    ip_run, olsen_run = ip["trajectory_run"], olsen["trajectory_run"]
    ip_seconds, olsen_seconds = float(ip_run["seconds"]), float(olsen_run["seconds"])
    if ip_run["correction_evaluations"] != len(rows) or olsen_run["correction_evaluations"] != len(rows):
        raise ValueError("correction report count mismatch")
    thresholds = []
    for tolerance in (1e-4, 1e-5, 1e-6):
        a = next((r for r in rows if abs(r["ip_corrected_energy"] - ip["reference_energy"]) <= tolerance), None)
        b = next((r for r in rows if abs(r["ip_olsen_corrected_energy"] - ip["reference_energy"]) <= tolerance), None)
        thresholds.append({
            "tolerance_hartree": tolerance,
            "ip_iteration": a["iteration"] if a else None,
            "ip_olsen_iteration": b["iteration"] if b else None,
            "ip_seconds": a["ip_wall_seconds"] if a else None,
            "ip_olsen_seconds": b["ip_olsen_wall_seconds"] if b else None,
            "wall_time_speedup": a["ip_wall_seconds"] / b["ip_olsen_wall_seconds"] if a and b else None,
        })
    summary = {
        "ip_seconds": ip_seconds,
        "ip_olsen_seconds": olsen_seconds,
        "olsen_extra_wall_seconds": olsen_seconds - ip_seconds,
        "olsen_relative_overhead": (olsen_seconds - ip_seconds) / ip_seconds,
        "overhead_definition": "(IP+Olsen total seconds - IP total seconds) / IP total seconds",
        "ip_report_seconds": ip_run["correction_seconds"],
        "ip_olsen_report_seconds": olsen_run["correction_seconds"],
        "report_timer_difference_seconds": olsen_run["correction_seconds"] - ip_run["correction_seconds"],
        "trajectory_points": len(rows),
        "max_variational_energy_difference_hartree": max(
            abs(r["ip_variational_energy"] - r["ip_olsen_variational_energy"]) for r in rows),
        "target_tolerances": thresholds,
        "timing_note": "One run per variant; cross-job differences also include node load and OpenMP variation. Report timers exclude per-update cache maintenance.",
    }
    return summary, rows


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("ip_json", type=Path)
    parser.add_argument("ip_olsen_json", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    ip = json.loads(args.ip_json.read_text())
    olsen = json.loads(args.ip_olsen_json.read_text())
    summary, rows = compare(ip, olsen)
    summary.update(ip_input=str(args.ip_json.resolve()), ip_olsen_input=str(args.ip_olsen_json.resolve()))
    args.output.mkdir(parents=True, exist_ok=True)
    (args.output / "06_ip_olsen_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    with (args.output / "06_ip_olsen_trajectory.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.2))
    iterations = [r["iteration"] for r in rows]
    axes[0].plot(iterations, [r["ip_wall_seconds"] for r in rows], label="IP only")
    axes[0].plot(iterations, [r["ip_olsen_wall_seconds"] for r in rows], label="IP + Olsen")
    axes[0].set(xlabel="iteration", ylabel="cumulative wall time (s)")
    axes[0].legend()
    axes[1].semilogx([r["olsen_stored_determinants"] for r in rows],
                    [100 * r["olsen_relative_overhead"] for r in rows])
    axes[1].set(xlabel="nonzero variational coefficients", ylabel="extra time relative to IP (%)")
    for ax in axes:
        ax.grid(True, alpha=0.25)
    fig.tight_layout()
    fig.savefig(args.output / "06_ip_olsen_timing.pdf")
    fig.savefig(args.output / "06_ip_olsen_timing.png", dpi=220)
    plt.close(fig)
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
