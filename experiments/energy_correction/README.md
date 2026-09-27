# Energy-correction paper experiments

This directory turns every TODO in the companion paper's
[`sections/04_numerical_results.tex`](../../../paper/sections/04_numerical_results.tex)
into a reproducible batch workflow. Submit from the CDFCI project root with:

```bash
sbatch experiments/energy_correction/run_all.sbatch
```

The production profile runs on the `bigMem3` partition for at most 18 hours,
requests 64 CPU cores and 96 GB, and does not request an exclusive node. OpenMP
threads are bound to physical cores. The default output is
`experiment_results/energy_correction/<job-id>/`. Use `small` or `large` as the
optional first argument to run only experiments 1--4 or experiment 5. For
the inexpensive experiments 1--4 without a whole-node allocation, use:

```bash
sbatch experiments/energy_correction/run_small.sbatch
```

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
   (C) local correction overhead versus stored wavefunction entries. Panel A
   uses the number of nonzero variational coefficients (`|x|_0`); panel C uses
   the number of stored residual/wavefunction entries (`|z|_0`), because those
   entries determine the correction scan cost. Local overhead is the correction
   evaluation time divided by the total wall time in the corresponding
   reporting interval. `05_trajectory.csv` contains every plotted value.
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
C2 and 220 N2 points while keeping repeated full-wavefunction correction scans
within the 18-hour allocation. One paired raw/corrected timing run is included.

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
