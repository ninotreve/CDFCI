#!/bin/bash
#SBATCH --partition=partition
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=128G
#SBATCH --time=18:00:00
#SBATCH --job-name=cdfci_h2o_corr
#SBATCH --output=cdfci_h2o_corr_%j.out
#SBATCH --error=cdfci_h2o_corr_%j.out

# One-job H2O/cc-pVDZ energy-correction experiment.
#
# The job runs, from fresh wavefunctions on the same allocation:
#   1. raw CDFCI and IP+Olsen (paired trajectory benchmark), then
#   2. IP-only (same solver settings), followed by analysis and plots.
#
# The default 4,000,000 block iterations is chosen to leave time for both
# correction variants, JSONL flushing, and post-processing within 18 hours.
# Override CDFCI_H2O_ITERATIONS after a short calibration run if desired.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="${CDFCI_PROJECT_DIR:-$(cd "${SCRIPT_DIR}/../.." && pwd)}"
BUILD_DIR="${CDFCI_EXPERIMENT_BUILD_DIR:-${PROJECT_DIR}/build_h2o_experiment}"
APP="${CDFCI_BENCHMARK_APP:-${BUILD_DIR}/bin/cdfci_energy_correction_benchmark_omp}"
FCIDUMP="${PROJECT_DIR}/data/h2o_ccpvdz_psi4.FCIDUMP"

ITERATIONS="${CDFCI_H2O_ITERATIONS:-4000000}"
INTERVAL="${CDFCI_CORRECTION_INTERVAL:-10000}"
Z_THRESHOLD="${CDFCI_H2O_Z_THRESHOLD:-0}"
MAX_WAVEFUNCTION_SIZE="${CDFCI_H2O_MAX_WAVEFUNCTION_SIZE:-848256040}"
NUM_COORDINATES="${CDFCI_NUM_COORDINATES:-64}"
TIMING_REPEATS="${CDFCI_TIMING_REPEATS:-1}"
REFERENCE_ENERGY="${CDFCI_H2O_REFERENCE_ENERGY:--76.2418600697}"

if [[ -n "${CDFCI_EXPERIMENT_OUTPUT:-}" ]]; then
    OUTPUT_ROOT="${CDFCI_EXPERIMENT_OUTPUT}"
elif [[ -n "${SCRATCH:-}" ]]; then
    OUTPUT_ROOT="${SCRATCH}/cdfci_h2o_energy_correction/${SLURM_JOB_ID:-manual}"
else
    OUTPUT_ROOT="${PROJECT_DIR}/experiment_results/energy_correction/h2o_ccpvdz/${SLURM_JOB_ID:-manual}"
fi
OUTPUT="${OUTPUT_ROOT}/h2o_ccpvdz"

if [[ -f /public3/soft/modules/module.sh ]]; then
    source /public3/soft/modules/module.sh
    module load gcc/12.2
    module load cmake/3.28.0
fi

export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-64}"
export OMP_DYNAMIC=false
export OMP_PROC_BIND=spread
export OMP_PLACES=cores
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export MPLBACKEND=Agg
export MPLCONFIGDIR="${SLURM_TMPDIR:-/tmp}/cdfci-h2o-mpl-${SLURM_JOB_ID:-manual}"
export CDFCI_BENCHMARK_VERBOSE=1
export REFERENCE_ENERGY_FOR_MANIFEST="${REFERENCE_ENERGY}"
export NUM_COORDINATES_FOR_MANIFEST="${NUM_COORDINATES}"
export ITERATIONS_FOR_MANIFEST="${ITERATIONS}"
export INTERVAL_FOR_MANIFEST="${INTERVAL}"
export Z_THRESHOLD_FOR_MANIFEST="${Z_THRESHOLD}"
export MAX_WF_FOR_MANIFEST="${MAX_WAVEFUNCTION_SIZE}"

[[ "${OMP_NUM_THREADS}" == "64" ]] || {
    echo "This experiment is planned for 64 OpenMP threads; got ${OMP_NUM_THREADS}." >&2
    exit 2
}
[[ -f "${FCIDUMP}" ]] || { echo "Missing FCIDUMP: ${FCIDUMP}" >&2; exit 1; }
[[ "${ITERATIONS}" =~ ^[0-9]+$ ]] && (( ITERATIONS > 0 )) || {
    echo "CDFCI_H2O_ITERATIONS must be a positive integer." >&2; exit 2;
}
[[ "${INTERVAL}" =~ ^[0-9]+$ ]] && (( INTERVAL > 0 )) || {
    echo "CDFCI_CORRECTION_INTERVAL must be a positive integer." >&2; exit 2;
}

mkdir -p "${OUTPUT}" "${MPLCONFIGDIR}"
cd "${PROJECT_DIR}"

echo "H2O/cc-pVDZ CDFCI energy-correction experiment"
echo "project=${PROJECT_DIR}"
echo "fcidump=${FCIDUMP}"
echo "output=${OUTPUT}"
echo "threads=${OMP_NUM_THREADS} coordinates=${NUM_COORDINATES}"
echo "iterations=${ITERATIONS} interval=${INTERVAL} z_threshold=${Z_THRESHOLD}"
echo "max_wavefunction_size=${MAX_WAVEFUNCTION_SIZE} reference=${REFERENCE_ENERGY}"
echo "started=$(date '+%F %T %z')"

# Configure and compile the exact checkout submitted with this job. Set
# CDFCI_EIGEN_SOURCE if compute nodes cannot download Eigen during configure.
if [[ -n "${CDFCI_BENCHMARK_APP:-}" ]]; then
    [[ -x "${APP}" ]] || { echo "Missing CDFCI_BENCHMARK_APP: ${APP}" >&2; exit 1; }
else
    CMAKE_ARGS=(
        -S "${PROJECT_DIR}"
        -B "${BUILD_DIR}"
        -DCMAKE_BUILD_TYPE=Release
        -DBUILD_TESTING=OFF
        -DCDFCI_ENABLE_OPENMP=ON
        -DCDFCI_BUILD_EXPERIMENTS=ON
    )
    if [[ -n "${CDFCI_EIGEN_SOURCE:-}" ]]; then
        CMAKE_ARGS+=("-DFETCHCONTENT_SOURCE_DIR_EIGEN=${CDFCI_EIGEN_SOURCE}")
    fi
    cmake "${CMAKE_ARGS[@]}"
    cmake --build "${BUILD_DIR}" \
        --target cdfci_energy_correction_benchmark_omp \
        --parallel "${CDFCI_BUILD_JOBS:-8}"
fi

[[ -x "${APP}" ]] || { echo "Benchmark executable not found: ${APP}" >&2; exit 1; }

REFERENCE_FILE="${OUTPUT}/reference.json"
python3 - "${REFERENCE_FILE}" <<'PY'
import json
import os
import sys
from pathlib import Path

Path(sys.argv[1]).write_text(json.dumps({
    "energy": float(os.environ["REFERENCE_ENERGY_FOR_MANIFEST"]),
    "source": "perturbation/CDFCI/regression_tests/h2o/ccpvdz_psi4/output",
    "method": "300000000-iteration CDFCI calculation",
    "fcidump": "data/h2o_ccpvdz_psi4.FCIDUMP",
    "note": "Long-run CDFCI reference; reference resolution must be checked against the trajectory tail."
}, indent=2) + "\n")
PY

python3 - "${OUTPUT}/manifest.json" "${APP}" "${FCIDUMP}" <<'PY'
import hashlib
import json
import os
import platform
import subprocess
import sys
from pathlib import Path

manifest_path, app, fcidump = sys.argv[1:]
source_paths = [
    "test/benchmark_energy_correction.cpp",
    "include/energy_correction.h",
    "include/ip_entry.h",
    "include/wavefunction.h",
    "include/solver.h",
    "experiments/energy_correction/run_h2o_ccpvdz.sbatch",
]
manifest = {
    "experiment": "h2o_ccpvdz_ip_ip_olsen",
    "git_commit": subprocess.check_output(["git", "rev-parse", "HEAD"], text=True).strip(),
    "job_id": os.environ.get("SLURM_JOB_ID"),
    "node": os.environ.get("SLURMD_NODENAME"),
    "partition": os.environ.get("SLURM_JOB_PARTITION"),
    "threads": int(os.environ["OMP_NUM_THREADS"]),
    "coordinates": int(os.environ["NUM_COORDINATES_FOR_MANIFEST"]),
    "iterations": int(os.environ["ITERATIONS_FOR_MANIFEST"]),
    "report_interval": int(os.environ["INTERVAL_FOR_MANIFEST"]),
    "z_threshold": float(os.environ["Z_THRESHOLD_FOR_MANIFEST"]),
    "max_wavefunction_size": int(os.environ["MAX_WF_FOR_MANIFEST"]),
    "reference_energy": float(os.environ["REFERENCE_ENERGY_FOR_MANIFEST"]),
    "execution_order": ["raw+ip+olsen", "ip-only"],
    "omp_proc_bind": os.environ["OMP_PROC_BIND"],
    "omp_places": os.environ["OMP_PLACES"],
    "fcidump": str(Path(fcidump).resolve()),
    "fcidump_sha256": hashlib.sha256(Path(fcidump).read_bytes()).hexdigest(),
    "benchmark_app": str(Path(app).resolve()),
    "benchmark_app_sha256": hashlib.sha256(Path(app).read_bytes()).hexdigest(),
    "source_sha256": {
        path: hashlib.sha256(Path(path).read_bytes()).hexdigest()
        for path in source_paths
    },
    "host": platform.platform(),
}
Path(manifest_path).write_text(json.dumps(manifest, indent=2) + "\n")
PY

PAIRED_JSON="${OUTPUT}/05_large_scale_with_timings.json"
IP_JSON="${OUTPUT}/05_ip_only.json"
ARGS=("${FCIDUMP}" "${ITERATIONS}" "${INTERVAL}" "${REFERENCE_ENERGY}")
# The final zero explicitly disables early stopping: both paired trajectories
# contain exactly the requested steps and the final partial report is retained.
SOLVER_ARGS=("${Z_THRESHOLD}" "${MAX_WAVEFUNCTION_SIZE}" "${NUM_COORDINATES}" "${TIMING_REPEATS}" "0" "${Z_THRESHOLD}" "0")

echo "Running raw + IP + Olsen paired benchmark"
srun --ntasks=1 --cpus-per-task="${OMP_NUM_THREADS}" --cpu-bind=cores \
    "${APP}" "${ARGS[@]}" "${PAIRED_JSON}" timed "${SOLVER_ARGS[@]}"

echo "Running IP-only benchmark"
srun --ntasks=1 --cpus-per-task="${OMP_NUM_THREADS}" --cpu-bind=cores \
    "${APP}" "${ARGS[@]}" "${IP_JSON}" ip-only "${SOLVER_ARGS[@]}"

python3 - "${PAIRED_JSON}" "${IP_JSON}" "${OUTPUT}" "${REFERENCE_ENERGY}" <<'PY'
import csv
import json
import math
import sys
from pathlib import Path

paired_path, ip_path, output_path, reference_text = sys.argv[1:]
paired = json.loads(Path(paired_path).read_text())
ip = json.loads(Path(ip_path).read_text())
output = Path(output_path)
reference = float(reference_text)
joint = paired["trajectory"]
ip_rows = ip["trajectory"]
if len(joint) != len(ip_rows):
    raise RuntimeError("IP and IP+Olsen trajectories have different report counts")

rows = []
for a, b in zip(joint, ip_rows):
    if a["iteration"] != b["iteration"]:
        raise RuntimeError("IP and IP+Olsen report iterations do not align")
    energies = {
        "variational": float(a["variational_energy"]),
        "ip": float(b["corrected_energy"]),
        "ip_olsen": float(a["corrected_energy"]),
    }
    row = {
        "iteration": a["iteration"],
        "variational_energy": energies["variational"],
        "ip_corrected_energy": energies["ip"],
        "ip_olsen_corrected_energy": energies["ip_olsen"],
        "variational_signed_error": energies["variational"] - reference,
        "ip_signed_error": energies["ip"] - reference,
        "ip_olsen_signed_error": energies["ip_olsen"] - reference,
        "variational_absolute_error": abs(energies["variational"] - reference),
        "ip_absolute_error": abs(energies["ip"] - reference),
        "ip_olsen_absolute_error": abs(energies["ip_olsen"] - reference),
        "ip_correction": float(b["external_correction"]),
        "olsen_internal_correction": float(a["internal_correction"]),
        "ip_wall_seconds": float(b["wall_seconds"]),
        "ip_olsen_wall_seconds": float(a["corrected_wall_seconds"]),
        "stored_determinants": a["stored_determinants"],
        "stored_wavefunction_entries": a["stored_wavefunction_entries"],
        "hamiltonian_columns": a["hamiltonian_columns"],
    }
    rows.append(row)

with (output / "h2o_signed_error_trajectory.csv").open("w", newline="") as stream:
    writer = csv.DictWriter(stream, fieldnames=list(rows[0]), lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)

def first_negative(key):
    return next((r["iteration"] for r in rows if r[key] < 0), None)

summary = {
    "reference_energy": reference,
    "trajectory_points": len(rows),
    "final_iteration": rows[-1]["iteration"],
    "final_signed_errors": {
        key: rows[-1][key]
        for key in ("variational_signed_error", "ip_signed_error", "ip_olsen_signed_error")
    },
    "first_below_reference_iteration": {
        "variational": first_negative("variational_signed_error"),
        "ip": first_negative("ip_signed_error"),
        "ip_olsen": first_negative("ip_olsen_signed_error"),
    },
    "ip_solver_seconds": ip["trajectory_run"]["seconds"],
    "ip_olsen_solver_seconds": paired["trajectory_run"]["seconds"],
    "ip_olsen_minus_ip_seconds": (
        paired["trajectory_run"]["seconds"] - ip["trajectory_run"]["seconds"]
    ),
}
(output / "h2o_signed_error_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
print(json.dumps(summary, indent=2))
PY

python3 experiments/energy_correction/summarize_large_scale.py \
    "${PAIRED_JSON}" "${OUTPUT}"
python3 experiments/energy_correction/compare_ip_olsen.py \
    "${IP_JSON}" "${PAIRED_JSON}" "${OUTPUT}"

python3 - "${OUTPUT}/h2o_signed_error_trajectory.csv" "${OUTPUT}" <<'PY'
import csv
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

source, output = sys.argv[1:]
with open(source, newline="") as stream:
    rows = list(csv.DictReader(stream))
iterations = [int(r["iteration"]) for r in rows]
fig, axes = plt.subplots(1, 2, figsize=(12, 4.5))
series = (
    ("variational_signed_error", "CDFCI", "tab:blue"),
    ("ip_signed_error", "CDFCI + IP", "tab:orange"),
    ("ip_olsen_signed_error", "CDFCI + IP + Olsen", "tab:green"),
)
for key, label, color in series:
    values = [float(r[key]) for r in rows]
    axes[0].plot(iterations, values, label=label, color=color, linewidth=1.2)
    axes[1].semilogy(iterations, [abs(value) for value in values], label=label,
                     color=color, linewidth=1.2)
for ax in axes:
    ax.axhline(0.0, color="black", linewidth=0.8)
    ax.grid(True, which="both", alpha=0.25)
axes[0].set(xlabel="block iteration", ylabel=r"signed error $E-E_{ref}$ (Ha)")
axes[0].set_yscale("symlog", linthresh=1e-10)
axes[1].set(xlabel="block iteration", ylabel=r"absolute error $|E-E_{ref}|$ (Ha)")
axes[0].legend()
axes[1].legend()
fig.tight_layout()
fig.savefig(Path(output) / "h2o_signed_error_curves.pdf")
fig.savefig(Path(output) / "h2o_signed_error_curves.png", dpi=220)
PY

echo "finished=$(date '+%F %T %z')"
echo "results=${OUTPUT}"
