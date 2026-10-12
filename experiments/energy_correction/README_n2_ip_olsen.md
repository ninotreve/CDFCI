# N2/cc-pVDZ IP+Olsen versus IP on amd_256

Submit `submit_n2_ip_olsen.slurm` from the repository root. It runs IP+Olsen,
then IP only, using fresh wavefunctions on the same allocated node. There are
64 OpenMP threads, 2200000 block steps, up to 64 coordinates per step, reports
every 10000 steps, `z_threshold=5e-7`, and an entry ceiling of 848256040.
These are PR #3's N2 settings. The script requests 128 GiB and 18 hours.

The benchmark takes the entry ceiling directly from the script; it does not
read an input JSON's `max_memory`. Slurm's `--mem` is the scheduler allocation.
The previous test with these parameters completed with 96 GiB allocated and a
recorded peak RSS of approximately 41.4 GiB. The new pending buffers hold only
changed internal determinant keys, with capacity retained between reports.

Compile separately after copying the updated sources, from the repository root:

```bash
source /public3/soft/modules/module.sh
module load gcc/12.2
module load cmake/3.28.0
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DBUILD_TESTING=OFF \
  -DCDFCI_ENABLE_OPENMP=ON -DCDFCI_BUILD_EXPERIMENTS=ON
cmake --build build --target cdfci_energy_correction_benchmark_omp --parallel 8
sbatch submit_n2_ip_olsen.slurm
```

Required input: `data/n2_ccpvdz_psi4.FCIDUMP`. GCC is loaded inside the job for
its runtime libraries. Python 3's standard library writes the comparison CSV;
NumPy and Matplotlib are not required.

The solver's usual energy rows are flushed to `n2_ip_olsen_<job-id>.out` during
each solve. Results are under
`experiment_results/energy_correction/<job-id>/n2_ccpvdz/`:

- `ip_olsen.json`: IP+Olsen trajectory and timing.
- `ip_only.json`: IP-only trajectory and timing.
- `runtime_summary.csv`: both total times, report correction times and final
  energies. `internal_b_refreshes` counts unique internal rows refreshed over
  all reports; `peak_pending_internal_rows` is the largest report queue.

The JSON for each invocation is saved after it finishes. If the second solve
is killed, the first JSON and the already flushed out rows remain available.
Thread scheduling can produce small trajectory differences between independent
runs. These two sequential runs provide one timing comparison; repetition on
the server is needed to establish a reliable speedup.

The pending strategy keeps the original Olsen formula and precision. It queues
each internal row at most once between reports, refreshes its final `b` through
the main table, then computes Olsen on the compact internal array with the
current `E0`. PT2 retains its existing incremental IP update.
