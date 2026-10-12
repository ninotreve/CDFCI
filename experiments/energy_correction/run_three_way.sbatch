#!/bin/bash
#SBATCH --partition=partition
#SBATCH --nodelist=bigMem5
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=96G
#SBATCH --time=18:00:00
#SBATCH --job-name=cdfci_three_way
#SBATCH --output=cdfci_three_way_%j.out
#SBATCH --error=cdfci_three_way_%j.out

set -euo pipefail

SYSTEM="${1:-}"
case "${SYSTEM}" in
    c2_ccpvdz)
        ITERATIONS="${CDFCI_C2_ITERATIONS:-800000}"
        Z_THRESHOLD="${CDFCI_C2_Z_THRESHOLD:-3e-8}"
        MAX_WAVEFUNCTION_SIZE="${CDFCI_C2_MAX_WAVEFUNCTION_SIZE:-1696512081}"
        REFERENCE_KEY=c2_ccpvdz_psi4
        ;;
    n2_ccpvdz)
        ITERATIONS="${CDFCI_N2_ITERATIONS:-2200000}"
        Z_THRESHOLD="${CDFCI_N2_Z_THRESHOLD:-5e-7}"
        MAX_WAVEFUNCTION_SIZE="${CDFCI_N2_MAX_WAVEFUNCTION_SIZE:-848256040}"
        REFERENCE_KEY=n2_ccpvdz_psi4
        ;;
    *) echo "Usage: sbatch $0 c2_ccpvdz|n2_ccpvdz" >&2; exit 2 ;;
esac

PROJECT_DIR="${CDFCI_PROJECT_DIR:-${SLURM_SUBMIT_DIR:-$PWD}}"
BUILD_DIR="${CDFCI_EXPERIMENT_BUILD_DIR:-${PROJECT_DIR}/build_experiments}"
APP="${CDFCI_BENCHMARK_APP:-${BUILD_DIR}/bin/cdfci_energy_correction_benchmark_omp}"
FCIDUMP="${PROJECT_DIR}/data/${SYSTEM}_psi4.FCIDUMP"
OUTPUT_ROOT="${CDFCI_EXPERIMENT_OUTPUT:-${PROJECT_DIR}/experiment_results/energy_correction/${SLURM_JOB_ID:-manual}}"
OUTPUT="${OUTPUT_ROOT}/${SYSTEM}"
INTERVAL="${CDFCI_CORRECTION_INTERVAL:-10000}"
NUM_COORDINATES="${CDFCI_NUM_COORDINATES:-64}"

if [[ -f /public3/soft/modules/module.sh ]]; then
    source /public3/soft/modules/module.sh
    module load gcc/12.2
fi

export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-64}"
export OMP_DYNAMIC=false
export OMP_PROC_BIND=spread
export OMP_PLACES=cores
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export MPLBACKEND=Agg
export MPLCONFIGDIR="${SLURM_TMPDIR:-/tmp}/cdfci-three-way-${SLURM_JOB_ID:-manual}"
# Every variant uses the same live report setting, including report I/O time.
export CDFCI_BENCHMARK_VERBOSE=1

cd "${PROJECT_DIR}"
[[ -x "${APP}" ]] || { echo "Missing benchmark executable: ${APP}" >&2; exit 1; }
[[ -f "${FCIDUMP}" ]] || { echo "Missing FCIDUMP: ${FCIDUMP}" >&2; exit 1; }
mkdir -p "${OUTPUT}" "${MPLCONFIGDIR}"

REFERENCE_ENERGY="$(python3 -c \
    'import json,sys; print(json.load(open(sys.argv[1]))[sys.argv[2]]["energy"])' \
    experiments/energy_correction/references.json "${REFERENCE_KEY}")"

python3 - "${OUTPUT_ROOT}" "${SYSTEM}" "${APP}" "${REFERENCE_KEY}" \
    "${ITERATIONS}" "${INTERVAL}" "${Z_THRESHOLD}" \
    "${MAX_WAVEFUNCTION_SIZE}" "${NUM_COORDINATES}" <<'PY'
import hashlib
import json
import os
import subprocess
import sys
from pathlib import Path

root, system, app, ref_key, iterations, interval, threshold, ceiling, coordinates = sys.argv[1:]
root = Path(root)
references = json.loads(Path('experiments/energy_correction/references.json').read_text())
(root / system / 'reference.json').write_text(json.dumps(references[ref_key], indent=2) + '\n')
manifest = {
    'mode': 'three-way-parallel-olsen',
    'benchmark_systems': [system],
    'git_commit': subprocess.check_output(['git', 'rev-parse', 'HEAD'], text=True).strip(),
    'job_id': os.environ.get('SLURM_JOB_ID'),
    'node': os.environ.get('SLURMD_NODENAME'),
    'partition': os.environ.get('SLURM_JOB_PARTITION'),
    'openmp_threads': int(os.environ['OMP_NUM_THREADS']),
    'omp_proc_bind': os.environ['OMP_PROC_BIND'],
    'omp_places': os.environ['OMP_PLACES'],
    'benchmark_app': str(Path(app).resolve()),
    'binary_sha256': hashlib.sha256(Path(app).read_bytes()).hexdigest(),
    'num_iterations': int(iterations),
    'report_interval': int(interval),
    'z_threshold': float(threshold),
    'max_wavefunction_size': int(ceiling),
    'num_coordinates': int(coordinates),
    'live_reports': True,
    'execution_order': ['raw', 'ip+olsen', 'ip'],
    'source_sha256': {
        path: hashlib.sha256(Path(path).read_bytes()).hexdigest()
        for path in ('test/benchmark_energy_correction.cpp', 'include/energy_correction.h',
                     'include/wavefunction.h', 'include/solver.h',
                     'experiments/energy_correction/run_three_way.sbatch')
    },
}
(root / 'manifest.json').write_text(json.dumps(manifest, indent=2) + '\n')
PY

echo "started=$(date '+%F %T %z')"
echo "system=${SYSTEM} partition=${SLURM_JOB_PARTITION:-unset} node=${SLURMD_NODENAME:-unset} threads=${OMP_NUM_THREADS} coordinates=${NUM_COORDINATES}"
echo "iterations=${ITERATIONS} interval=${INTERVAL} z_threshold=${Z_THRESHOLD} max_wavefunction_size=${MAX_WAVEFUNCTION_SIZE} reference_energy=${REFERENCE_ENERGY}"
echo "app=${APP} output=${OUTPUT}"

PAIRED_JSON="${OUTPUT}/05_large_scale_with_timings.json"
IP_JSON="${OUTPUT}/05_ip_only.json"
ARGS=("${FCIDUMP}" "${ITERATIONS}" "${INTERVAL}" "${REFERENCE_ENERGY}")
SOLVER_ARGS=("${Z_THRESHOLD}" "${MAX_WAVEFUNCTION_SIZE}" "${NUM_COORDINATES}" 1)

# Same-node raw and IP+Olsen, from separate fresh wavefunctions.
# This mode also runs the benchmark's standard 200-step warm-up.
srun --ntasks=1 --cpus-per-task="${OMP_NUM_THREADS}" --cpu-bind=cores \
    "${APP}" "${ARGS[@]}" "${PAIRED_JSON}" timed "${SOLVER_ARGS[@]}"

# Same-node IP-only, again from a fresh wavefunction with the same warm-up.
srun --ntasks=1 --cpus-per-task="${OMP_NUM_THREADS}" --cpu-bind=cores \
    "${APP}" "${ARGS[@]}" "${IP_JSON}" ip-only "${SOLVER_ARGS[@]}"

python3 - "${PAIRED_JSON}" "${IP_JSON}" "${OUTPUT}" <<'PY'
import csv
import json
import sys
from pathlib import Path

paired = json.loads(Path(sys.argv[1]).read_text())
ip = json.loads(Path(sys.argv[2]).read_text())
output = Path(sys.argv[3])
rows = [
    {'mode': label, 'solver_seconds': run['seconds'],
     'variational_energy': run['energy'], 'correction_seconds': run['correction_seconds'],
     'correction_evaluations': run['correction_evaluations']}
    for label, run in [('raw', paired['baseline_trajectory_run']),
                       ('ip', ip['trajectory_run']), ('ip+olsen', paired['trajectory_run'])]
]
assert rows[0]['correction_evaluations'] == 0
assert ip['options']['energy_correction']['enabled']
assert not ip['options']['energy_correction']['olsen_enabled']
assert paired['options']['energy_correction']['olsen_enabled']
assert all(r['internal_correction'] == 0 for r in ip['trajectory'])
with (output / 'runtime_summary.csv').open('w', newline='') as handle:
    writer = csv.DictWriter(handle, fieldnames=list(rows[0]), lineterminator='\n')
    writer.writeheader()
    writer.writerows(rows)
summary = {
    'runs': rows,
    'ip_overhead_vs_raw': rows[1]['solver_seconds'] / rows[0]['solver_seconds'] - 1,
    'olsen_overhead_vs_ip': rows[2]['solver_seconds'] / rows[1]['solver_seconds'] - 1,
    'total_overhead_vs_raw': rows[2]['solver_seconds'] / rows[0]['solver_seconds'] - 1,
    'timing_note': 'Same allocation, one solve per variant; report I/O included. Fixed order raw, IP+Olsen, IP.',
}
(output / 'runtime_summary.json').write_text(json.dumps(summary, indent=2) + '\n')
print(json.dumps(summary, indent=2))
PY

if python3 -c 'import numpy, matplotlib' >/dev/null 2>&1; then
    srun --ntasks=1 --cpus-per-task=1 \
        python3 -u experiments/energy_correction/summarize_large_scale.py \
        "${PAIRED_JSON}" "${OUTPUT}"
    srun --ntasks=1 --cpus-per-task=1 \
        python3 -u experiments/energy_correction/compare_ip_olsen.py \
        "${IP_JSON}" "${PAIRED_JSON}" "${OUTPUT}"
else
    echo "Plotting skipped: NumPy/Matplotlib unavailable; raw JSON and runtime summaries saved."
fi

echo "finished=$(date '+%F %T %z') results=${OUTPUT_ROOT}"
