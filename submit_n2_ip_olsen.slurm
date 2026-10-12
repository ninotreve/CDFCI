#!/bin/bash
#SBATCH -p amd_256
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=128G
#SBATCH -J n2_ip_olsen
#SBATCH -t 18:00:00
#SBATCH --output=n2_ip_olsen_%j.out
#SBATCH --error=n2_ip_olsen_%j.out

set -euo pipefail

source /public3/soft/modules/module.sh
module load gcc/12.2

# Compile separately and submit from the project root.
PROJECT_DIR="${CDFCI_PROJECT_DIR:-${SLURM_SUBMIT_DIR:-$PWD}}"
BUILD_DIR="${CDFCI_EXPERIMENT_BUILD_DIR:-${PROJECT_DIR}/build}"
APP="${BUILD_DIR}/bin/cdfci_energy_correction_benchmark_omp"
FCIDUMP="${PROJECT_DIR}/data/n2_ccpvdz_psi4.FCIDUMP"
OUTPUT="${CDFCI_EXPERIMENT_OUTPUT:-${PROJECT_DIR}/experiment_results/energy_correction/${SLURM_JOB_ID:-manual}}/n2_ccpvdz"

# Reproduce PR #3: 2,200,000 block steps, up to 64 coordinates per step.
ITERATIONS="${CDFCI_N2_ITERATIONS:-2200000}"
INTERVAL="${CDFCI_CORRECTION_INTERVAL:-10000}"
Z_THRESHOLD="${CDFCI_N2_Z_THRESHOLD:-5e-7}"
MAX_WAVEFUNCTION_SIZE="${CDFCI_N2_MAX_WAVEFUNCTION_SIZE:-848256040}"
NUM_COORDINATES="${CDFCI_NUM_COORDINATES:-64}"
REFERENCE_ENERGY=-109.2821730115

export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-64}"
export OMP_DYNAMIC=false
export OMP_PROC_BIND=spread
export OMP_PLACES=cores
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export CDFCI_BENCHMARK_VERBOSE=1

cd "${PROJECT_DIR}"
[[ -x "${APP}" ]] || { echo "Benchmark executable not found: ${APP}" >&2; exit 1; }
[[ -f "${FCIDUMP}" ]] || { echo "FCIDUMP not found: ${FCIDUMP}" >&2; exit 1; }
command -v python3 >/dev/null || { echo "Python 3 is required for the timing summary" >&2; exit 1; }
mkdir -p "${OUTPUT}"

echo "N2/cc-pVDZ threads=${OMP_NUM_THREADS} coordinates=${NUM_COORDINATES} iterations=${ITERATIONS} report_interval=${INTERVAL} z_threshold=${Z_THRESHOLD}"
echo "results=${OUTPUT}"

# Each invocation starts with a fresh wavefunction on the same allocation.
for MODE in olsen-only ip-only; do
    if [[ "${MODE}" == olsen-only ]]; then
        RESULT_FILE="${OUTPUT}/ip_olsen.json"
    else
        RESULT_FILE="${OUTPUT}/ip_only.json"
    fi
    srun --ntasks=1 --cpus-per-task="${OMP_NUM_THREADS}" --cpu-bind=cores \
        "${APP}" "${FCIDUMP}" "${ITERATIONS}" "${INTERVAL}" "${REFERENCE_ENERGY}" \
        "${RESULT_FILE}" "${MODE}" \
        "${Z_THRESHOLD}" "${MAX_WAVEFUNCTION_SIZE}" "${NUM_COORDINATES}" 1
done

# Standard library only. The ordinary solver out columns stay unchanged.
python3 - "${OUTPUT}" <<'PY'
import csv
import json
import pathlib
import sys

directory = pathlib.Path(sys.argv[1])
rows = []
for mode, filename in (("IP+Olsen", "ip_olsen.json"), ("IP", "ip_only.json")):
    with (directory / filename).open() as stream:
        result = json.load(stream)
    run = result["trajectory_run"]
    final = result["trajectory"][-1]
    rows.append({
        "mode": mode,
        "seconds": run["seconds"],
        "correction_seconds": run["correction_seconds"],
        "final_iteration": final["iteration"],
        "E0": final["variational_energy"],
        "PT2": final["external_correction"],
        "Olsen": final["internal_correction"] if mode == "IP+Olsen" else "",
        "Ecorr": final["corrected_energy"],
        "V": final["stored_determinants"],
        "Z": final["stored_wavefunction_entries"],
        "internal_b_refreshes": run.get("internal_b_refreshes", ""),
        "peak_pending_internal_rows": run.get("peak_pending_internal_rows", ""),
    })
with (directory / "runtime_summary.csv").open("w", newline="") as stream:
    writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
    writer.writeheader()
    writer.writerows(rows)
full, ip = rows
difference = full["seconds"] - ip["seconds"]
print("IP+Olsen: {:.2f} s; IP: {:.2f} s".format(full["seconds"], ip["seconds"]))
print("IP+Olsen minus IP: {:.2f} s ({:+.2f}%)".format(
    difference, 100 * difference / ip["seconds"]))
print("Final E0 difference: {:.12e} Ha".format(full["E0"] - ip["E0"]))
PY

echo "Completed: ${OUTPUT}"
