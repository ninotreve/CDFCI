# Energy-correction paper experiments

Completed C2/N2 runs, figures, timing ablations, and paper-ready findings are
indexed in [RESULTS.md](RESULTS.md), including limitations of the "almost free"
claim. The latest parallel-Olsen three-way results are documented in
[RESULTS_PARALLEL.md](RESULTS_PARALLEL.md): C2 job `1007783` and N2 job
`1007784`, with full-solve Olsen overheads of 5.14% and 23.21% over IP-only.
Only the six completed production result directories are versioned;
local tests, other runs, build files, and new FCIDUMP files remain ignored.

For a new same-node raw / IP-only / parallel IP+Olsen comparison on this
cluster's `partition` partition, use `run_three_way.sbatch`. It keeps the
previous 64-thread, 64-coordinate production settings, requests 96 GiB and
18 hours on bigMem5 without exclusive allocation, and runs all three variants
within one job per molecule. Build the latest OpenMP benchmark first:

```bash
cmake --build build_experiments --target cdfci_energy_correction_benchmark_omp --parallel 4
sbatch --job-name=cdfci_c2_three_way experiments/energy_correction/run_three_way.sbatch c2_ccpvdz
sbatch --job-name=cdfci_n2_three_way experiments/energy_correction/run_three_way.sbatch n2_ccpvdz
```

The execution order is raw, IP+Olsen, then IP-only, using a fresh wavefunction
for every solve. All variants use live reporting. Outputs include the original
JSON trajectories, `runtime_summary.csv` / JSON with all three solve times,
the existing raw/corrected and IP/Olsen comparison plots, and an early-written
manifest containing binary/source SHA256 hashes and exact solver parameters.
`CDFCI_BENCHMARK_APP` can select a frozen binary for reproducible queued jobs.
Olsen's incremental overhead is `(T_IP+Olsen - T_IP) / T_IP`; each variant has
one full measurement, so node activity and run order can still affect timings.

This directory turns every TODO in the companion paper's
[`sections/04_numerical_results.tex`](../../../paper/sections/04_numerical_results.tex)
into a reproducible batch workflow. Submit from the CDFCI project root with:

```bash
sbatch experiments/energy_correction/run_all.sbatch
```

The production profile runs on the `bigMem5` node for at most 18 hours,
requests 64 CPU cores and 96 GB, and does not request an exclusive node. OpenMP
threads are bound to physical cores. The default output is
`experiment_results/energy_correction/<job-id>/`. Use `small` or `large` as the
optional first argument to run only experiments 1--4 or experiment 5. The
`baseline2` mode runs IP/PT2 only (`enabled=true`, `olsen_enabled=false`) once
per system, with the same production parameters and warm-up as `large`. For
the inexpensive experiments 1--4 without a whole-node allocation, use:

```bash
sbatch experiments/energy_correction/run_small.sbatch
```

For production-scale experiment 5, submit C2 and N2 as independent jobs so
each system receives its own 18-hour allocation:

```bash
sbatch --job-name=cdfci_c2_corr --export=ALL,CDFCI_BENCHMARK_SYSTEMS=c2_ccpvdz \
  experiments/energy_correction/run_all.sbatch large
sbatch --job-name=cdfci_n2_corr --export=ALL,CDFCI_BENCHMARK_SYSTEMS=n2_ccpvdz \
  experiments/energy_correction/run_all.sbatch large
```

To measure the marginal cost of Olsen using an existing IP+Olsen JSON, submit
each system separately. Replace the JSON path below with that system's previous
`05_large_scale_with_timings.json`:

```bash
sbatch --job-name=cdfci_c2_ip \
  --export=ALL,CDFCI_BENCHMARK_SYSTEMS=c2_ccpvdz,CDFCI_IP_OLSEN_REFERENCE_JSON=/absolute/path/to/c2/05_large_scale_with_timings.json \
  experiments/energy_correction/run_all.sbatch baseline2
sbatch --job-name=cdfci_n2_ip \
  --export=ALL,CDFCI_BENCHMARK_SYSTEMS=n2_ccpvdz,CDFCI_IP_OLSEN_REFERENCE_JSON=/absolute/path/to/n2/05_large_scale_with_timings.json \
  experiments/energy_correction/run_all.sbatch baseline2
```

`05_ip_only.json` retains the IP-only trajectory and complete solve time.
`06_ip_olsen_summary.json`, CSV, PDF, and PNG compare matched iterations and
compute `(IP+Olsen time - IP time) / IP time`, including Olsen's internal-cache
maintenance. Configuration and Hamiltonian work must match. The report timer
difference is also exported, but excludes per-update maintenance. Cross-job
timings are single measurements and can be affected by node load; they should
not alone establish that Olsen is "almost free".


## Experiment map

1. `01_synthetic_orders.*`: the fixed 3-by-3 model `H=D+epsilon B`, fitted on
   the eight smallest couplings. It reports the predicted orders 2, 3, and 4.
2. `02_davidson_trajectory.*`: Davidson on an explicitly specified 400-by-400
   half-filled, open six-site Hubbard CI Hamiltonian, including energy errors,
   residual norms, and the postprocessed/raw error ratio at every iterate.
3. `03_poor_initial_vector.*`: the same Hamiltonian from a seeded poor random
   vector, retaining every iterate and the first ratio below one.
4. `04_hubbard_failure_sweep.*`: a `U/t` sweep at the deliberately early first
   Ritz iterate, including error-ratio degradation and overshoot flags.
5. `05_*`: OpenMP block-CDFCI trajectories for **C2/cc-pVDZ and N2/cc-pVDZ**,
   using the repository FCIDUMPs and the production settings recorded in
   `/home/yjzhang/codes/cdfci/logs/c2n2-bm3.out`: 64 threads, 64 coordinates,
   `z_threshold=3e-8` for C2 and `5e-7` for N2. C2 runs 800,000 iterations with
   a 1,696,512,081-entry ceiling; N2 runs 2,200,000 iterations with an
   848,256,040-entry ceiling. Accuracy is measured against the repository's
   much longer CDFCI regression energies (`-75.7319603747` for C2 and
   `-109.2821730115` for N2), without rerunning a reference calculation. The
   finite-iteration endpoints from `c2n2-bm3.out` remain recorded in
   `references.json`, but are not suitable error references because they retain
   the truncation error that the correction is intended to estimate.

   The generated `05_correction_advantage.pdf` and PNG contain the requested
   panels: (A) error versus stored determinants, (B) error versus wall time, and
   (C) cumulative end-to-end correction overhead versus stored wavefunction
   entries. Panel A uses the number of nonzero variational coefficients
   (`|x|_0`); panel C uses the number of stored residual/wavefunction entries
   (`|z|_0`), because the streaming PT2 work and diagonal cache scale with the
   stored wavefunction. The plotted overhead is derived from paired trajectories
   as `(corrected wall time - raw wall time) / corrected wall time`, so it
   includes both report-time Olsen evaluation and per-update incremental PT2
   cost. The direct `correction_seconds` timer is retained as a diagnostic only.
   `05_trajectory.csv` contains every plotted value.
   Because the available long-run references do not resolve the sub-`1e-7` Ha
   regime, panels A and B omit all error points below `1e-7` Ha and state this
   reference-resolution floor explicitly. Panel C remains reference-independent.

The experiments require Python, NumPy, and Matplotlib. CMake needs Eigen;
if it is not installed and compute nodes cannot download it, point at an
unpacked Eigen source tree:

```bash
sbatch --export=ALL,CDFCI_EIGEN_SOURCE=/path/to/eigen-3.4.0 \
  experiments/energy_correction/run_all.sbatch
```

The correction/reporting interval defaults to 10,000 iterations. This gives 80
C2 and 220 N2 points. The optimized implementation maintains PT2 incrementally
and recomputes Olsen over the compact internal array at report points. One
paired raw/corrected timing run is included.

Useful global overrides are `CDFCI_EXPERIMENT_OUTPUT`,
`CDFCI_BENCHMARK_SYSTEMS`, `CDFCI_BENCHMARK_ITERATIONS`,
`CDFCI_CORRECTION_INTERVAL`, `CDFCI_Z_THRESHOLD`,
`CDFCI_MAX_WAVEFUNCTION_SIZE`, `CDFCI_NUM_COORDINATES`,
`CDFCI_TIMING_REPEATS`, and `CDFCI_OMP_THREADS`. Per-system iteration and
storage overrides use the `CDFCI_C2_*` and `CDFCI_N2_*` prefixes shown in
`run_all.sbatch`. The iteration count should be an integer multiple of the
correction interval so plots have uniform sampling.

A custom FCIDUMP can still be run by setting both `CDFCI_BENCHMARK_FCIDUMP` and
`CDFCI_REFERENCE_ENERGY`. The workflow never attempts infeasible exact
diagonalization automatically.
