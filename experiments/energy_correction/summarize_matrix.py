#!/usr/bin/env python3
"""Combine experiment-5 system summaries into manuscript-ready tables."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("output_root", type=Path)
    parser.add_argument("systems", nargs="+")
    args = parser.parse_args()

    systems = []
    tolerance_rows = []
    for name in args.systems:
        path = args.output_root / name / "05_large_scale_summary.json"
        summary = json.loads(path.read_text(encoding="utf-8"))
        system_row = {
            "system": name,
            "spatial_orbitals": summary["spatial_orbitals"],
            "spin_orbitals": summary["spin_orbitals"],
            "electrons": summary["electrons"],
            "fixed_ms_determinant_dimension": summary["fixed_ms_determinant_dimension"],
            "reference_energy": summary["reference_energy"],
            "reference_iterations": summary["reference_iterations"],
            "reference_tail_energy_change": summary["reference_tail_energy_change"],
            "improved_point_fraction": summary["improved_point_fraction"],
            "overshoot_points": summary["overshoot_points"],
            "invalid_correction_points": summary["invalid_correction_points"],
            "direct_relative_runtime_overhead": summary["direct_relative_runtime_overhead"],
            "median_paired_runtime_overhead": summary["median_relative_runtime_overhead"],
        }
        systems.append(system_row)
        for threshold in summary["target_tolerances"]:
            tolerance_rows.append({"system": name, **threshold})

    matrix = {"systems": systems, "target_tolerances": tolerance_rows}
    (args.output_root / "05_benchmark_matrix.json").write_text(
        json.dumps(matrix, indent=2) + "\n", encoding="utf-8")
    for filename, rows in (("05_systems.csv", systems),
                           ("05_target_tolerances_all_systems.csv", tolerance_rows)):
        if not rows:
            continue
        with (args.output_root / filename).open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
            writer.writeheader()
            writer.writerows(rows)
    print(json.dumps(matrix, indent=2))


if __name__ == "__main__":
    main()
