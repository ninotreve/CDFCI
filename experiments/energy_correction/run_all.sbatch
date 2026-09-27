#!/bin/bash
#SBATCH -w bigMem3
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=96G
#SBATCH -J cdfci_pe
#SBATCH -t 18:00:00
#SBATCH --output=cdfci_pe_%j.out
#SBATCH --error=cdfci_pe_%j.out

set -euo pipefail

# Project directory:
#   1. Use CDFCI_PROJECT_DIR if explicitly provided.
#   2. Under Slurm, use the directory from which sbatch was submitted.
#   3. Otherwise, use the current working directory.
PROJECT_DIR="${CDFCI_PROJECT_DIR:-${SLURM_SUBMIT_DIR:-$PWD}}"

MODE="${1:-all}"
case "${MODE}" in
    all|small|large) ;;
    *) echo "Usage: sbatch $0 [all|small|large]" >&2; exit 2 ;;
esac

if [[ -f /public3/soft/modules/module.sh ]]; then
    source /public3/soft/modules/module.sh
    module load gcc/12.2
    module load cmake/3.28.0
fi

THREADS="${SLURM_CPUS_PER_TASK:-64}"
export OMP_NUM_THREADS="${CDFCI_OMP_THREADS:-${THREADS}}"
export OMP_DYNAMIC=false
export OMP_PROC_BIND=spread
export OMP_PLACES=cores
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export PYSCF_NUM_THREADS=1

export MPLCONFIGDIR="${SLURM_TMPDIR:-/tmp}/cdfci-matplotlib-${SLURM_JOB_ID:-manual}"
mkdir -p "${MPLCONFIGDIR}"

OUTPUT_ROOT="${CDFCI_EXPERIMENT_OUTPUT:-${PROJECT_DIR}/experiment_results/energy_correction/${SLURM_JOB_ID:-manual}}"
mkdir -p "${OUTPUT_ROOT}"

MANIFEST_SYSTEMS=()

echo "CDFCI energy-correction experiment suite"
echo "mode=${MODE}"
echo "project=${PROJECT_DIR}"
echo "output=${OUTPUT_ROOT}"
echo "openmp_threads=${OMP_NUM_THREADS}"
echo "slurm_submit_dir=${SLURM_SUBMIT_DIR:-unset}"
echo "started=$(date '+%F %T %z')"

cd "${PROJECT_DIR}"

python3 -c 'import numpy, matplotlib; print("numpy", numpy.__version__, "matplotlib", matplotlib.__version__)'

if [[ "${MODE}" = all || "${MODE}" = small ]]; then
    srun --ntasks=1 --cpus-per-task=1 \
        python3 experiments/energy_correction/run_experiments.py \
        --output "${OUTPUT_ROOT}/small"
fi

if [[ "${MODE}" = all || "${MODE}" = large ]]; then
    BUILD_DIR="${CDFCI_EXPERIMENT_BUILD_DIR:-${PROJECT_DIR}/build_experiments}"

    cmake_args=(
        -S "${PROJECT_DIR}"
        -B "${BUILD_DIR}"
        -DCMAKE_BUILD_TYPE=Release
        -DBUILD_TESTING=OFF
        -DCDFCI_ENABLE_OPENMP=ON
        -DCDFCI_BUILD_EXPERIMENTS=ON
    )

    if [[ -n "${CDFCI_EIGEN_SOURCE:-}" ]]; then
        cmake_args+=("-DFETCHCONTENT_SOURCE_DIR_EIGEN=${CDFCI_EIGEN_SOURCE}")
    fi

    cmake "${cmake_args[@]}"

    BUILD_JOBS="${CDFCI_BUILD_JOBS:-8}"

    cmake --build "${BUILD_DIR}" \
        --target cdfci_energy_correction_benchmark_omp \
        --parallel "${BUILD_JOBS}"

    if [[ -n "${CDFCI_BENCHMARK_FCIDUMP:-}" ]]; then
        if [[ -z "${CDFCI_REFERENCE_ENERGY:-}" ]]; then
            echo "A custom FCIDUMP requires CDFCI_REFERENCE_ENERGY." >&2
            exit 2
        fi

        BENCHMARK_SYSTEMS=(custom)
    else
        read -r -a BENCHMARK_SYSTEMS <<< \
            "${CDFCI_BENCHMARK_SYSTEMS:-c2_ccpvdz n2_ccpvdz}"
    fi

    MANIFEST_SYSTEMS=("${BENCHMARK_SYSTEMS[@]}")

    for SYSTEM in "${BENCHMARK_SYSTEMS[@]}"; do
        case "${SYSTEM}" in
            c2_ccpvdz)
                FCIDUMP="${PROJECT_DIR}/data/c2_ccpvdz_psi4.FCIDUMP"
                ITERATIONS="${CDFCI_C2_ITERATIONS:-${CDFCI_BENCHMARK_ITERATIONS:-800000}}"
                Z_THRESHOLD="${CDFCI_C2_Z_THRESHOLD:-${CDFCI_Z_THRESHOLD:-3e-8}}"
                MAX_WAVEFUNCTION_SIZE="${CDFCI_C2_MAX_WAVEFUNCTION_SIZE:-${CDFCI_MAX_WAVEFUNCTION_SIZE:-1696512081}}"
                REFERENCE_KEY=c2_ccpvdz_psi4
                ;;

            n2_ccpvdz)
                FCIDUMP="${PROJECT_DIR}/data/n2_ccpvdz_psi4.FCIDUMP"
                ITERATIONS="${CDFCI_N2_ITERATIONS:-${CDFCI_BENCHMARK_ITERATIONS:-2200000}}"
                Z_THRESHOLD="${CDFCI_N2_Z_THRESHOLD:-${CDFCI_Z_THRESHOLD:-5e-7}}"
                MAX_WAVEFUNCTION_SIZE="${CDFCI_N2_MAX_WAVEFUNCTION_SIZE:-${CDFCI_MAX_WAVEFUNCTION_SIZE:-848256040}}"
                REFERENCE_KEY=n2_ccpvdz_psi4
                ;;

            custom)
                FCIDUMP="${CDFCI_BENCHMARK_FCIDUMP}"
                ITERATIONS="${CDFCI_BENCHMARK_ITERATIONS:-200000}"
                Z_THRESHOLD="${CDFCI_Z_THRESHOLD:-5e-5}"
                MAX_WAVEFUNCTION_SIZE="${CDFCI_MAX_WAVEFUNCTION_SIZE:-848256040}"
                REFERENCE_ENERGY="${CDFCI_REFERENCE_ENERGY}"
                REFERENCE_KEY=""
                ;;

            *)
                echo "Unknown benchmark system: ${SYSTEM}" >&2
                exit 2
                ;;
        esac

        if [[ ! -f "${FCIDUMP}" ]]; then
            echo "FCIDUMP not found: ${FCIDUMP}" >&2
            exit 1
        fi

        INTERVAL="${CDFCI_CORRECTION_INTERVAL:-10000}"
        NUM_COORDINATES="${CDFCI_NUM_COORDINATES:-64}"
        TIMING_REPEATS="${CDFCI_TIMING_REPEATS:-1}"

        SYSTEM_OUTPUT="${OUTPUT_ROOT}/${SYSTEM}"
        mkdir -p "${SYSTEM_OUTPUT}"

        if [[ -n "${REFERENCE_KEY}" ]]; then
            REFERENCE_ENERGY="$(python3 -c \
                'import json,sys; print(json.load(open(sys.argv[1]))[sys.argv[2]]["energy"])' \
                "${PROJECT_DIR}/experiments/energy_correction/references.json" \
                "${REFERENCE_KEY}")"
            python3 - "${PROJECT_DIR}/experiments/energy_correction/references.json" \
                "${REFERENCE_KEY}" "${SYSTEM_OUTPUT}/reference.json" <<'PY'
import json
import sys
from pathlib import Path

references = json.loads(Path(sys.argv[1]).read_text(encoding="utf-8"))
Path(sys.argv[3]).write_text(
    json.dumps(references[sys.argv[2]], indent=2) + "\n", encoding="utf-8")
PY
        fi

        TIMING_JSON="${SYSTEM_OUTPUT}/05_large_scale_with_timings.json"

        echo "system=${SYSTEM} fcidump=${FCIDUMP} iterations=${ITERATIONS} reference_energy=${REFERENCE_ENERGY}"

        srun --ntasks=1 \
            --cpus-per-task="${OMP_NUM_THREADS}" \
            --cpu-bind=cores \
            "${BUILD_DIR}/bin/cdfci_energy_correction_benchmark_omp" \
            "${FCIDUMP}" \
            "${ITERATIONS}" \
            "${INTERVAL}" \
            "${REFERENCE_ENERGY}" \
            "${TIMING_JSON}" \
            timed \
            "${Z_THRESHOLD}" \
            "${MAX_WAVEFUNCTION_SIZE}" \
            "${NUM_COORDINATES}" \
            "${TIMING_REPEATS}"

        srun --ntasks=1 --cpus-per-task=1 \
            python3 experiments/energy_correction/summarize_large_scale.py \
            "${TIMING_JSON}" \
            "${SYSTEM_OUTPUT}"
    done

    srun --ntasks=1 --cpus-per-task=1 \
        python3 experiments/energy_correction/summarize_matrix.py \
        "${OUTPUT_ROOT}" \
        "${BENCHMARK_SYSTEMS[@]}"
fi

python3 - "${OUTPUT_ROOT}" "${MODE}" "${MANIFEST_SYSTEMS[*]}" <<'PY'
import json
import os
import platform
import subprocess
import sys
from pathlib import Path

root, mode, systems = Path(sys.argv[1]), sys.argv[2], sys.argv[3].split()

try:
    commit = subprocess.check_output(
        ["git", "rev-parse", "HEAD"],
        text=True
    ).strip()
except Exception:
    commit = None

manifest = {
    "mode": mode,
    "benchmark_systems": systems,
    "git_commit": commit,
    "job_id": os.environ.get("SLURM_JOB_ID"),
    "node": os.environ.get("SLURMD_NODENAME"),
    "python": platform.python_version(),
    "openmp_threads": int(os.environ.get("OMP_NUM_THREADS", "1")),
    "slurm_cpus_per_task": os.environ.get("SLURM_CPUS_PER_TASK"),
    "slurm_submit_dir": os.environ.get("SLURM_SUBMIT_DIR"),
    "output_root": str(root),
}

(root / "manifest.json").write_text(
    json.dumps(manifest, indent=2) + "\n"
)
PY

echo "finished=$(date '+%F %T %z')"
echo "results=${OUTPUT_ROOT}"
