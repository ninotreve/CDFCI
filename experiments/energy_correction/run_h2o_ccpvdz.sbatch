#!/bin/bash
#SBATCH --partition=partition
#SBATCH --nodelist=bigMem3
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=128G
#SBATCH --time=18:00:00
#SBATCH --signal=B:USR1@600
#SBATCH --job-name=cdfci_h2o_5m
#SBATCH --output=cdfci_h2o_5m_%j.out
#SBATCH --error=cdfci_h2o_5m_%j.out

# Compile and seal a unique source snapshot BEFORE submission; see README_h2o.md.
# Exactly one 5M-block IP+Olsen solve; three energies share the same state.
# No calibration, warmups, independent baselines, repetitions or step reduction.
set -euo pipefail
case "${1:-}" in
    ''|--check-only) ;;
    *) echo "Usage: $0 [--check-only]" >&2; exit 2 ;;
esac
if [[ -n "${CDFCI_PROJECT_DIR:-}" ]]; then
    PROJECT_DIR="${CDFCI_PROJECT_DIR}"
elif [[ -n "${SLURM_SUBMIT_DIR:-}" ]]; then
    PROJECT_DIR="${SLURM_SUBMIT_DIR}"
else
    SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
    PROJECT_DIR="$(cd "${SCRIPT_DIR}/../.." && pwd)"
fi
[[ -f "${PROJECT_DIR}/CMakeLists.txt" && -f "${PROJECT_DIR}/test/benchmark_energy_correction.cpp" ]] || {
    echo "Invalid CDFCI project directory: ${PROJECT_DIR}. Submit from the repository root or set CDFCI_PROJECT_DIR." >&2
    exit 2
}
PROJECT_DIR="$(cd "${PROJECT_DIR}" && pwd)"
ITERATIONS=5000000
INTERVAL="${CDFCI_CORRECTION_INTERVAL:-10000}"
Z_THRESHOLD="${CDFCI_H2O_Z_THRESHOLD:-0}"
MAX_WAVEFUNCTION_SIZE="${CDFCI_H2O_MAX_WAVEFUNCTION_SIZE:-848256040}"
NUM_COORDINATES="${CDFCI_NUM_COORDINATES:-64}"
REFERENCE_ENERGY=-76.2418600697
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-64}"
export OMP_DYNAMIC=false OMP_PROC_BIND=spread OMP_PLACES=cores
export OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
export CDFCI_BENCHMARK_VERBOSE=1 MPLBACKEND=Agg
[[ "${OMP_NUM_THREADS}" == 64 && "${NUM_COORDINATES}" == 64 ]] || {
    echo 'This experiment requires 64 threads and 64 coordinates.' >&2; exit 2;
}
[[ -f "${PROJECT_DIR}/data/h2o_ccpvdz_psi4.FCIDUMP" ]] || {
    echo "Missing FCIDUMP: ${PROJECT_DIR}/data/h2o_ccpvdz_psi4.FCIDUMP" >&2; exit 1;
}
[[ "${ITERATIONS}" =~ ^[0-9]+$ && "${INTERVAL}" =~ ^[0-9]+$ ]] &&
    (( ITERATIONS > 0 && INTERVAL > 0 )) || { echo 'Invalid iteration settings' >&2; exit 2; }
if [[ "${1:-}" == --check-only ]]; then
    echo "project=${PROJECT_DIR}"
    echo "fcidump=${PROJECT_DIR}/data/h2o_ccpvdz_psi4.FCIDUMP"
    echo "iterations=${ITERATIONS} interval=${INTERVAL} z_threshold=${Z_THRESHOLD} (fixed length; no automatic reduction)"
    echo "threads=${OMP_NUM_THREADS} coordinates=${NUM_COORDINATES} max_wavefunction_size=${MAX_WAVEFUNCTION_SIZE} reference=${REFERENCE_ENERGY}"
    echo 'solves=1 mode=IP+Olsen calibration=0 timing_pairs=0 node=bigMem3'
    echo 'Path and parameter checks passed; no build or solver was started.'
    exit 0
fi
FROZEN="${CDFCI_FROZEN_RUN:-${PROJECT_DIR}/build_h2o_experiment/frozen-h2o-single-20261019}"
[[ -f "${FROZEN}/source_manifest.json" && -x "${FROZEN}/build/bin/cdfci_energy_correction_benchmark_omp" ]]
[[ -f "${FROZEN}/src/experiments/energy_correction/h2o_single_run.py" && -f "${FROZEN}/src/reference/h2o_ccpvdz_new.out" ]] || {
    echo 'A sealed single-run snapshot including the certified reference is required; see README_h2o.md.' >&2; exit 2;
}
if [[ -f /public3/soft/modules/module.sh ]]; then
    source /public3/soft/modules/module.sh
    module load gcc/12.2
fi
OUTPUT_ROOT="${CDFCI_EXPERIMENT_OUTPUT:-${PROJECT_DIR}/experiment_results/energy_correction/h2o_ccpvdz/${SLURM_JOB_ID:?Submit with sbatch}}"
export MPLCONFIGDIR="${SLURM_TMPDIR:-/tmp}/cdfci-h2o-mpl-${SLURM_JOB_ID}"
mkdir -p "${MPLCONFIGDIR}"
echo "started=$(date '+%F %T %z') frozen=${FROZEN} output=${OUTPUT_ROOT} node=${SLURMD_NODENAME:-unset}"
# exec lets Slurm's pre-timeout USR1 reach the runner, which saves partial
# analysis and terminates only its own srun process group before the hard limit.
exec python3 "${FROZEN}/src/experiments/energy_correction/h2o_single_run.py" run \
    --frozen "${FROZEN}" --output "${OUTPUT_ROOT}" \
    --iterations "${ITERATIONS}" --interval "${INTERVAL}" \
    --z-threshold "${Z_THRESHOLD}" \
    --max-wavefunction-size "${MAX_WAVEFUNCTION_SIZE}" \
    --reference-energy "${REFERENCE_ENERGY}"
