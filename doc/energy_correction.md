# Streaming energy correction during CDFCI

The optional `solver.cdfci.energy_correction` postprocessor implements the
projected diagonal/Olsen energy estimator in *Low-Cost Perturbative Energy
Correction for Unconverged Configuration Interaction Iterates*, Algorithm 1,
Eqs. (15)-(19). It uses stored coefficients `x`, stored `z`, and cached
Hamiltonian diagonal elements. It never applies H to a new vector, changes the
CDFCI coefficients, or changes coordinate selection or stopping criteria.

## Enable

Add this object under `solver.cdfci` in an existing input:

```json
"energy_correction": {
    "enabled": true,
    "olsen_enabled": true,
    "store_history": false
}
```

The default is disabled. Corrections are reported at the solver's
`report_interval`; an early-convergence or final snapshot is also taken when
needed, without duplicating a report-point snapshot. If a run stops before
its first report, correction columns remain zero. Single-state CDFCI supports
both single-coordinate and multicoordinate updates, with serial and OpenMP
storage. Multiple excited states are rejected when correction is enabled.

The main wavefunction stores `c+d`, `b`, and a cached `H_ii` for every stored
row. The diagonal is evaluated when a row first enters the table and reused
when its external PT2 amplitude changes. PT2 keeps the incremental IP update
strategy: untouched rows retain their earlier denominator. The ownership bit
stored with `H_ii` distinguishes `c` from `d` without a second hash lookup.
When Olsen is enabled, a separate compact internal-only array mirrors `c`,
`b`, and `H_ii`; Olsen is recomputed from the current Rayleigh energy at report
points without scanning the main table. The internal determinant index is
used only to maintain that compact array.

The main progress table contains `PT2`, `Olsen`, and `Ecorr` columns. `PT2` is
`external_correction`, `Olsen` is `internal_correction`, and `Ecorr` is the
variational energy plus both corrections. When the estimator is disabled,
correction columns are zero. The top-level `report_interval` in the input
controls both progress output and the full Olsen recomputation schedule. PT2
tracking starts from the initial wavefunction and continues between reports.

Run `demo_input_energy_correction.json` from the `examples` directory. Its input
path is relative to that directory. There is no top-level `perturbation` section
in this example, so the old end-of-run Hx reconstruction is not requested.

## Relation to the predecessor and definition of the spaces

The predecessor, *Perturbative Coordinate Descent Full Configuration Interaction*
(DOI 10.1021/acs.jctc.6c00103), uses the compressed space already stored by CDFCI.
Here `internal_correction` applies the new projected estimator within supp(x).
`external_correction` is the predecessor's scalar PT2 sum on **stored** entries
with x_i = 0 and z_i != 0. The direct snapshot evaluator uses the current
Rayleigh energy for every denominator. The IP-CDFCI streaming path updates a
stored external coefficient when its z entry changes, using that update's
Rayleigh energy; unchanged entries keep their previous denominator until they
change again. This makes streaming PT2 an incremental approximation to a full
current-energy snapshot.

The supported modes are:

* Olsen plus PT2: set `enabled: true` and `olsen_enabled: true` (the default).
* PT2 only: set `enabled: true` and `olsen_enabled: false`.
* No correction: set `enabled: false`. The underlying CDFCI iteration is
  unchanged.

PT2 is always included when energy correction is enabled because Olsen must not
be reported without PT2. The corrected energy is E + `internal_correction` +
`external_correction`.

The decomposition is algebraically exact for the supplied x, z and diagonal:
for x_i = 0, q_i = -z_i / ||x||, so the external entries affect only S_qq.
Their contribution is sum z_i^2 / ((E - H_ii) * (x'x)). Do not add the old PT2
result again to `corrected_energy`. The old top-level `perturbation` calculation
remains a separate, explicitly requested operation, with its own reconstruction
cost and potentially different external space.

## Cost and precision

The report-time Olsen evaluator scans the compact internal array twice and
does not traverse the much larger external hash table. Its cost and additional
memory are O(number of internal determinants). The main table uses one more
double per stored row for the diagonal cache. The external PT2 result comes
from the incremental IP accumulator; each changed external row reuses its
diagonal after the first evaluation. PT2 denominators for untouched rows can lag behind the
current energy by design. `report_interval` controls how often Olsen is
recomputed; it does not refresh unchanged PT2 denominators. Denominators are used
directly without gap checks or status reports; zero or small gaps can therefore
produce nonfinite or inaccurate corrections.

Accumulator arithmetic uses the same `QUAD_PRECISION` type as CDFCI (GCC
`__float128`, or the existing long-double fallback). The implementation handles
unnormalised, signed x and the multicoordinate scaling convention. Eliminating
the pivot allows E = H_pp, including the single-determinant HF case. The
calculation does not audit denominator gaps or reject nonfinite results; the
reported values are the direct arithmetic result.

Under normal consistent initialisation, CDFCI's exact selected-coordinate
recalculation and subsequent updates preserve z_i = (Hx)_i on supp(x), up to
roundoff, even with compression. Outside supp(x), stored z may be inaccurate and
other residual entries may be absent. Therefore the external result under
compression is an approximation on the stored space, not full-space EN-PT2.
Checkpoint data must correspond to the same Hamiltonian and energy origin.
The three-value wavefunction record is incompatible with checkpoints written
by the preceding two-value build; start a new run when upgrading that build.
`compressed_z` reports whether the configured solver compression threshold is
positive; it is not an independent certification of residual accuracy.

The paper's third-order accuracy statement requires its nearly diagonal
assumptions and accurate residual. It is not a universal improvement guarantee
for compressed residuals or strongly off-diagonal Hamiltonians. Corrected energies
can overshoot below the exact ground-state energy and are not variational bounds.

## Results and Python

`Result.energy` retains the original solver energy. New fields are:

* `energy_correction`: last snapshot, including iteration, variational energy,
  internal/external correction, corrected energy, residual norms, diagonal
  evaluation count, and elapsed evaluator time. Legacy validity/status fields
  remain available in the API and are not used to audit denominator gaps.
* `energy_correction_evaluations`, `energy_correction_seconds`: totals for the
  retained solver attempt (automatic threshold restarts reset them).
* `energy_correction_history`: snapshots, only when `store_history: true`.
  This optional history uses O(number of snapshots) memory, not O(determinants).

The Python `CDFCIResult` exposes these same fields. History collection works with
`verbose: 0`. Without history only the final snapshot and totals are retained.

## Validation and benchmarking

`cdfci_energy_correction_test` checks the streaming formula against a direct
projected solve, internal/external decomposition, scale and shift invariance,
the single-determinant HF case, and the near-diagonal orders 2/3/4.
It also compares original and corrected CDFCI trajectories on the repository's
H2O fixture, including compressed internal z and multicoordinate updates.
CTest registers serial and available OpenMP versions.

`test/benchmark_energy_correction.cpp` underlies both serial and OpenMP
benchmark targets. Compile with the project's usual include paths and
Eigen/quadmath, then pass FCIDUMP, iteration count, correction interval,
reference energy (or `auto`) and output JSON path. It supports up to 128 spin
orbitals, records an accuracy trajectory, and reuses that baseline/enabled pair
as the first timing pair. Loading FCIDUMP is outside the timings; solver
allocation is inside. `auto` first runs the requested longer reference budget
and records its tail energy change. A supplied reference must match the actual
input Hamiltonian, not merely the molecule name.
