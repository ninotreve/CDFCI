# Streaming energy correction during CDFCI

The optional `solver.cdfci.energy_correction` postprocessor implements the
projected diagonal/Olsen energy estimator in *Low-Cost Perturbative Energy
Correction for Unconverged Configuration Interaction Iterates*, Algorithm 1,
Eqs. (15)-(19). It uses only stored coefficients `x`, stored `z`, and on-demand
Hamiltonian diagonal elements. It never applies H to a new vector, changes the
CDFCI coefficients, or changes coordinate selection or stopping criteria.

## Enable

Add this object under `solver.cdfci` in an existing input:

```json
"energy_correction": {
    "enabled": true,
    "scope": "stored",
    "interval": 1000,
    "store_history": false,
    "gap_tolerance": 1e-12
}
```

The default is disabled. `interval: 0` uses `report_interval`. A snapshot is also
taken at the final iteration, including early convergence; duplicate final
snapshots are avoided. Single-state CDFCI supports both single-coordinate and
multicoordinate updates, with serial and OpenMP storage. Multiple excited states
are rejected when correction is enabled.

The main progress table always contains `PT2` and `E_corrected` columns so that
cluster logs have a stable format. `PT2` is `external_correction`, while
`E_corrected` is E + `internal_correction` + `external_correction`. When the
estimator is disabled, both columns are zero. For one correction at every
printed progress point, set `energy_correction.interval` equal to
`report_interval` (for example, set both to 1000).

Run `demo_input_energy_correction.json` from the `examples` directory. Its input
path is relative to that directory. There is no top-level `perturbation` section
in this example, so the old end-of-run Hx reconstruction is not requested.

## Relation to the predecessor and definition of the spaces

The predecessor, *Perturbative Coordinate Descent Full Configuration Interaction*
(DOI 10.1021/acs.jctc.6c00103), uses the compressed space already stored by CDFCI.
Here `internal_correction` applies the new projected estimator within supp(x).
`external_correction` is the predecessor's scalar PT2 sum on **stored** entries
with x_i = 0 and z_i != 0, evaluated with the current energy at each snapshot.
It does not use IP-CDFCI's stale per-coordinate perturbative coefficients.

* `scope: "internal"`: only relax the current nonzero coefficient space. External
  entries are skipped. The corrected energy is E + delta_internal.
* `scope: "stored"` (default): include the existing external entries as well.
  The corrected energy is E + delta_internal + delta_external.

The decomposition is algebraically exact for the supplied x, z and diagonal:
for x_i = 0, q_i = -z_i / ||x||, so the external entries affect only S_qq.
Their contribution is sum z_i^2 / ((E - H_ii) * (x'x)). Do not add the old PT2
result again to `corrected_energy`. The old top-level `perturbation` calculation
remains a separate, explicitly requested operation, with its own reconstruction
cost and potentially different external space.

## Cost and precision

The evaluator stores only a pivot and scalar accumulators. There are no new
per-determinant arrays, hash entries, wavefunctions, or diagonal caches. A normal
evaluation uses one coefficient scan and one diagonal scan, with one diagonal
evaluation per visited determinant. A problematic nonpivot internal denominator
can cause one additional diagonal scan after changing the pivot.

Cost is O(number of stored entries), plus diagonal evaluation cost. This is not
O(1) time per CDFCI iteration: scanning a large hash table every iteration can be
expensive. Choose a correction interval appropriate to the reporting/accuracy
needs and inspect the measured `energy_correction_seconds`. Diagonal elements
are deliberately recomputed rather than cached to preserve memory usage.

Accumulator arithmetic uses the same `QUAD_PRECISION` type as CDFCI (GCC
`__float128`, or the existing long-double fallback). The implementation handles
unnormalised, signed x and the multicoordinate scaling convention. Eliminating
the pivot allows E = H_pp, including the single-determinant HF case. Small
nonpivot gaps or a nonpositive/nearly singular projected denominator produce an
invalid status and NaN corrected energy; they are not silently clipped.

Under normal consistent initialisation, CDFCI's exact selected-coordinate
recalculation and subsequent updates preserve z_i = (Hx)_i on supp(x), up to
roundoff, even with compression. Outside supp(x), stored z may be inaccurate and
other residual entries may be absent. Therefore the external result under
compression is an approximation on the stored space, not full-space EN-PT2.
Checkpoint data must correspond to the same Hamiltonian and energy origin.
`compressed_z` reports whether the configured solver compression threshold is
positive; it is not an independent certification of residual accuracy.

The paper's third-order accuracy statement requires its nearly diagonal
assumptions and accurate residual. It is not a universal improvement guarantee
for compressed residuals or strongly off-diagonal Hamiltonians. Corrected energies
can overshoot below the exact ground-state energy and are not variational bounds.

## Results and Python

`Result.energy` retains the original solver energy. New fields are:

* `energy_correction`: last snapshot, including iteration, variational energy,
  internal/external correction, corrected energy, validity/status, residual norms,
  diagonal evaluation count, and elapsed evaluator time.
* `energy_correction_evaluations`, `energy_correction_seconds`: totals for the
  retained solver attempt (automatic threshold restarts reset them).
* `energy_correction_history`: snapshots, only when `store_history: true`.
  This optional history uses O(number of snapshots) memory, not O(determinants).

The Python `CDFCIResult` exposes these same fields. History collection works with
`verbose: 0`. Without history only the final snapshot and totals are retained.
For an internal-only calculation, `stored_residual_norm` covers the chosen
internal scope; no omitted external residual norm is inferred.

## Validation and benchmarking

`cdfci_energy_correction_test` checks the streaming formula against a direct
projected solve, internal/external decomposition, scale and shift invariance,
pivot singularities, invalid denominators, and the near-diagonal orders 2/3/4.
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
