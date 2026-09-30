#!/bin/bash
#SBATCH -p amd_256
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH -J c2_ccpvdz_4c
#SBATCH -t 12:00:00
#SBATCH --output=c2_ccpvdz_4c_%j.out
#SBATCH --error=c2_ccpvdz_4c_%j.out

set -euo pipefail

source /public3/soft/modules/module.sh
module load gcc/12.2

export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK}"
export OMP_DYNAMIC=false
export OMP_PROC_BIND=close
export OMP_PLACES=cores

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="${CDFCI_PROJECT_DIR:-${SLURM_SUBMIT_DIR:-${SCRIPT_DIR}}}"
APP="${PROJECT_DIR}/build/bin/cdfci_omp"
INPUT_ARG="${1:-examples/c2_ccpvdz_100k_pt_4c.json}"

if [[ "${INPUT_ARG}" = /* ]]; then
    INPUT="${INPUT_ARG}"
else
    INPUT="${PROJECT_DIR}/${INPUT_ARG}"
fi

if [[ ! -f "${INPUT}" ]]; then
    echo "Error: input JSON not found: ${INPUT}" >&2
    exit 1
fi

if [[ ! -x "${APP}" ]]; then
    echo "Error: compiled OpenMP executable not found: ${APP}" >&2
    exit 1
fi

# Resolve the FCIDUMP path relative to the JSON directory.
cd "$(dirname "${INPUT}")"
srun --ntasks=1 --cpus-per-task="${SLURM_CPUS_PER_TASK}" \
     "${APP}" "$(basename "${INPUT}")" 2>&1 \
    | sed -u -e '/^\[Energy correction\]/d' \
          -e '/^Energy correction evaluations:/d'
