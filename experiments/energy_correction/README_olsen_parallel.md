# C2 parallel Olsen test on amd_256

Copy these files over the matching paths in the server checkout:

- `include/energy_correction.h`
- `include/wavefunction.h`
- `test/benchmark_energy_correction.cpp`
- `submit_c2_olsen.slurm`

Compile from the project root (not inside the submission script):

```bash
source /public3/soft/modules/module.sh
module load gcc/12.2
module load cmake/3.28.0
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DBUILD_TESTING=OFF -DCDFCI_ENABLE_OPENMP=ON -DCDFCI_BUILD_EXPERIMENTS=ON
cmake --build build --target cdfci_energy_correction_benchmark_omp --parallel 8
sbatch submit_c2_olsen.slurm
```

The script runs only IP/PT2+Olsen with the previous C2 parameters: 800,000
block steps, up to 64 coordinates per step, a report every 10,000 steps,
`z_threshold=3e-8`, 1,696,512,081 maximum wavefunction entries, 64 CPUs, and
96 GB of Slurm memory. The entry limit and Slurm memory request are different
settings, just as in the previous benchmark.

Change `#SBATCH --cpus-per-task=64` to change the CPU count. Both CDFCI and
Olsen use that count through `OMP_NUM_THREADS`. For example:

```bash
sbatch --cpus-per-task=4 submit_c2_olsen.slurm
```

Changing the CPU count does not change the 64-coordinate block size.
The report interval and other solver parameters are near the top of the
script and can also be overridden by the existing `CDFCI_*` environment
variables.

Olsen still reads the compact internal array and selects the pivot serially.
The expensive sum uses a static OpenMP partition and private 128-bit sums,
followed by a reduction. Reports with fewer than 4,096 rows stay serial.
The formula, precision, PT2 update, and wavefunction update are unchanged.
The different summation order can change the last rounding bits.

The internal determinant-to-row index uses the existing flat robinhood map.
Entries are inserted before parallel b updates. An internal b update now queues
the determinant only on its first change since the previous report, using the
main entry's metadata to deduplicate across threads. At report time, queued
rows are refreshed by key in parallel before the Olsen sum. Repeated updates
do not query the internal index. The main wavefunction table and external PT2
diagonal cache keep their existing size and numeric values.

Normal energy rows are flushed to `c2_olsen_omp_JOBID.out` during the solve.
The final solver time is printed there. On successful completion, the JSON
is saved under `experiment_results/energy_correction/JOBID/c2_ccpvdz/olsen_parallel.json`.
No NumPy, Matplotlib, plotting step, or baseline solve is required.
