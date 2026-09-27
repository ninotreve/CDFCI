# Core-code audit (2026-09-24)

Scope: solver control flow, wavefunction/container storage, energy correction,
program drivers, Python result propagation, build/test wiring, and the numerical
experiment helpers. The audit used warning-enabled GCC builds plus the serial
and OpenMP CTest suites.

## Findings fixed in this change

1. **High — shifted energies were wrong at `verbose: 0`.**
   `CDFCISolver::output_final` returned before restoring the Hamiltonian shift,
   so programmatic and benchmark results could use a different energy origin.
   Shift restoration now precedes the verbosity-only output guard.
2. **High — OptOrbFCI resume used uninitialized state.**
   `optimal_orbitals.init_iteration != 0` entered an empty history-loading
   branch and then consumed uninitialized `energy`, `zero_rdm`, and RDM arrays.
   Resume now fails explicitly until history loading is implemented.
3. **Medium — requested iterations were silently truncated.**
   The main loop ran `num_iterations / report_interval` complete batches and
   discarded a remainder. It now executes a final partial report batch.
4. **Medium — quiet runs lost all trajectory metadata.**
   Energy, determinant-count, Hamiltonian-column-work, and wall-time histories
   were appended only after the `verbose == 0` return or not recorded at all.
   They are now recorded independently of printing, cleared between restarted
   attempts, and returned with `report_interval`.
5. **Medium — the simple Python facade discarded correction results.**
   `CDFCI(...).run()` copied only `Result.energy`. It now returns the complete
   `Result`, including correction history and timing totals.
6. **Low — the documented benchmark was not a CMake target.**
   `CDFCI_BUILD_EXPERIMENTS=ON` now exposes
   `cdfci_energy_correction_benchmark`, allowing the Slurm workflow to build it
   without hand-written compiler commands.

## Open findings

1. **High — concurrent-container iterators escape their lock.**
   `ContainerCuckoo::{begin,end,find}` creates a `locked_table`, obtains an
   iterator, unlocks it, and returns the iterator. Any use can race or dereference
   an iterator whose guard lifetime has ended. Internal read-only code mostly
   uses the safe `loop` wrapper, but `WaveFunctionBase::calculate_xx_xz` and the
   public iterator API remain unsafe for the concurrent container. The durable
   fix is an RAII locked-range API, not another unlocked iterator wrapper.
2. **Medium — binary checkpoints are not portable or self-describing.**
   Raw `size_t`, determinant, floating-point, and array representations are
   written without a magic number, schema version, endianness, type widths, or
   checksum. A mismatched build can silently misread a checkpoint before a
   stream error is detected.
3. **Medium — exact trajectory equivalence is not promised by OpenMP runs.**
   Parallel first insertion and hash-table ordering can retain different small
   compressed `z` entries. Existing tests correctly check invariants within one
   parallel trajectory, but publication timing/accuracy comparisons should use
   the serial benchmark or report run-to-run variability.
4. **Low — warning debt obscures higher-value compiler diagnostics.**
   The warning-enabled build emits many signed/unsigned comparisons, deprecated
   implicit copy-assignment warnings for `Determinant`, unused parameters, and a
   backslash in a line comment. None failed the current tests, but a staged
   warning cleanup is advisable before enabling `-Werror`.

## Remaining test coverage gaps

- No regression test exercises a positive initial energy and verifies shift
  restoration with `verbose: 0`.
- No test covers OptOrbFCI checkpoint/resume because the feature is not
  implemented.
- Checkpoint corruption, cross-compiler compatibility, and concurrent iterator
  misuse are untested.

The energy-correction regression now covers a `1001`-iteration quiet run with a
`100`-iteration report interval, including the final partial batch and all five
trajectory-history vectors.
