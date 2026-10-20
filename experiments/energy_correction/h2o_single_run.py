#!/usr/bin/env python3
"""One fixed-length H2O IP+Olsen solve, with three energies on the same state."""
from __future__ import annotations

import argparse
import importlib
import json
import math
import os
from pathlib import Path
import shutil
import signal
import subprocess
import time

from h2o_experiment import csv_file, digest, journal, prepare, read, save, seal, verify

STOP = False
ENERGIES = ('variational', 'ip', 'ip_olsen')


def prepare_single(args):
    reference = args.reference_file.resolve()
    if not reference.is_file():
        raise FileNotFoundError(reference)
    prepare(args.root, args.frozen, args.eigen)
    destination = args.frozen / 'src/reference/h2o_ccpvdz_new.out'
    destination.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(reference, destination)
    manifest = read(args.frozen / 'source_manifest.json')
    manifest['files']['src/reference/h2o_ccpvdz_new.out'] = digest(destination)
    manifest['reference_original_path'] = str(reference)
    save(args.frozen / 'source_manifest.json', manifest)


def analyze(output, plots=False):
    config = read(output / 'manifest.json')
    events = journal(output / 'runs/joint.json.jsonl')
    rows = [r for r in events if r.get('record_type') == 'report']
    starts = [r for r in events if r.get('record_type') == 'run_start']
    ends = [r for r in events if r.get('record_type') == 'run_end']
    if starts:
        if len(starts) != 1 or len(ends) > 1:
            raise ValueError('Expected one joint solve')
        options = starts[0]['options']
        flags = options['energy_correction']
        if not flags['enabled'] or not flags['olsen_enabled']:
            raise ValueError('Both IP and Olsen must be enabled')
        if options['num_iterations'] != config['iterations'] or options['stopping_dx_threshold'] != 0:
            raise ValueError('Fixed iteration count / disabled early stopping required')
        if options['z_threshold'] != config['z_threshold']:
            raise ValueError('Residual threshold differs from manifest')
    trajectory = []
    previous_iteration, previous_time = 0, -1
    ref = config['reference_energy']
    for r in rows:
        if r['phase'] != 'ip+olsen' or r['iteration'] <= previous_iteration or r['wall_seconds'] <= previous_time:
            raise ValueError('Invalid phase or non-increasing trajectory')
        previous_iteration, previous_time = r['iteration'], r['wall_seconds']
        energy = r['variational_energy']
        external, internal = r['external_correction'], r['internal_correction']
        values = {'variational': energy, 'ip': energy + external, 'ip_olsen': r['corrected_energy']}
        if not all(math.isfinite(x) for x in [*values.values(), external, internal, r['wall_seconds']]):
            raise ValueError('Non-finite report')
        if abs(values['ip_olsen'] - (energy + external + internal)) > 1e-12:
            raise ValueError('Inconsistent energy decomposition')
        row = {'iteration': r['iteration'], 'wall_seconds_joint': r['wall_seconds'],
               'nonzero_coefficients': r['stored_determinants'],
               'stored_wavefunction_entries': r['stored_wavefunction_entries'],
               'hamiltonian_columns': r['hamiltonian_columns'],
               'external_correction': external, 'olsen_correction': internal,
               'direct_report_seconds': r['correction_seconds']}
        for name, value in values.items():
            row[name + '_energy'] = value
            row[name + '_signed_error'] = value - ref
        trajectory.append(row)
    complete = bool(ends and rows and rows[-1]['iteration'] == config['iterations']
                    and ends[0]['summary']['iterations'] == config['iterations'])
    targets = []
    for name in ENERGIES:
        for tol in (1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8, 1e-9, 1e-10):
            hit = next((r for r in trajectory if abs(r[name + '_signed_error']) <= tol), None)
            targets.append({'energy_estimate': name, 'tolerance_hartree': tol,
                            'first_recorded_iteration': hit['iteration'] if hit else None,
                            'wall_seconds_joint': hit['wall_seconds_joint'] if hit else None})
    summary = {'experiment': 'single_joint', 'planned_iterations': config['iterations'],
               'last_iteration': rows[-1]['iteration'] if rows else 0,
               'reports': len(rows), 'complete': complete,
               'reference_policy': config['reference_policy'],
               'latest': trajectory[-1] if trajectory else None,
               'below_reference_reports': {name: sum(r[name + '_signed_error'] < 0 for r in trajectory)
                                           for name in ENERGIES},
               'solve_summary': ends[0]['summary'] if ends else None,
               'note': 'All three estimates use the same wavefunction and joint-run clock. No standalone method runtimes or performance comparisons. First recorded hits are sampled, not interpolated crossings.'}
    csv_file(output / 'three_way_trajectory.csv', trajectory)
    csv_file(output / 'target_tolerances.csv', targets)
    save(output / 'summary.json', summary)
    if plots and trajectory:
        make_plots(output, trajectory)
    return summary


def make_plots(output, rows):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    colors = ('#3c566e', '#dc872a', '#14846b')
    labels = ('Variational E', 'E + IP', 'E + IP + Olsen')
    for filename, xkeys, xlabels in (
            ('signed_errors', ('iteration', 'wall_seconds_joint'), ('Block iteration', 'Joint-run wall time (s)')),
            ('errors_vs_work', ('hamiltonian_columns', 'nonzero_coefficients'),
             ('Hamiltonian columns', 'Nonzero variational coefficients'))):
        fig, axes = plt.subplots(1, 2, figsize=(12, 4.6))
        errors = [r[name + '_signed_error'] for r in rows for name in ENERGIES]
        for ax, key, xlabel in zip(axes, xkeys, xlabels):
            for name, label, color in zip(ENERGIES, labels, colors):
                ax.plot([r[key] for r in rows], [r[name + '_signed_error'] for r in rows],
                        '.-', label=label, color=color, markersize=2)
            ax.axhline(0, color='black', linewidth=.7)
            ax.set_yscale('symlog', linthresh=5e-11)
            ax.set_ylim(min(-5e-11, min(errors)*1.15), max(5e-11, max(errors)*1.15))
            ax.set(xlabel=xlabel, ylabel='Signed E - E_ref (Ha)')
            ax.legend(fontsize=9)
            ax.grid(True, which='both', alpha=.2)
        fig.suptitle(f'H2O/cc-pVDZ: one IP+Olsen solve, same wavefunction; last block {rows[-1]["iteration"]:,}')
        fig.tight_layout()
        for ext in ('png', 'pdf'):
            fig.savefig(output / (filename + '.' + ext), dpi=180)
        plt.close(fig)


def stop_signal(signum, frame):
    global STOP
    STOP = True


def run(args):
    global STOP
    STOP = False
    frozen, output = args.frozen.resolve(), args.output.resolve()
    source = verify(frozen, True)
    for dependency in ('numpy', 'matplotlib'):
        importlib.import_module(dependency)
    reference = frozen / 'src/reference/h2o_ccpvdz_new.out'
    if not reference.is_file():
        raise FileNotFoundError('Prepare a new single-run snapshot including the certified reference')
    output.mkdir(parents=True, exist_ok=True)
    if (output / 'manifest.json').exists():
        raise RuntimeError('Refusing to overwrite an existing experiment')
    config = {k: str(v) if isinstance(v, Path) else v for k, v in vars(args).items()}
    config.update(experiment='single_joint', source=source,
                  job_id=os.environ.get('SLURM_JOB_ID'), node=os.environ.get('SLURMD_NODENAME'),
                  omp_threads=os.environ.get('OMP_NUM_THREADS'),
                  batch_process_affinity_cpu_ids=sorted(os.sched_getaffinity(0)),
                  omp_places=os.environ.get('OMP_PLACES'), omp_proc_bind=os.environ.get('OMP_PROC_BIND'),
                  reference_policy={'reference_energy': args.reference_energy,
                                    'source': str(reference), 'sha256': digest(reference),
                                    'original_path': source.get('reference_original_path'),
                                    'certification_provenance': 'User confirms the historical reference is strictly certified.',
                                    'rounding_half_unit_hartree': 5e-11,
                                    'note': 'No arbitrary 1e-7 exclusion region. Plot symlog linear scale uses printing resolution; it is not a total-energy certification bound.'},
                  timing_note='64 threads; physical-core count is not asserted. All three energies share one joint-run clock.')
    save(output / 'manifest.json', config)
    save(output / 'plan.json', {'iterations': args.iterations, 'number_of_solves': 1,
                              'calibration_solves': 0, 'warmup_solves': 0, 'timing_pairs': 0,
                              'automatic_iteration_reduction': False})
    target = output / 'runs/joint.json'
    target.parent.mkdir()
    app = frozen / 'build/bin/cdfci_energy_correction_benchmark_omp'
    cmd = ['srun', '--ntasks=1', '--cpus-per-task=64', '--cpu-bind=verbose,cores', str(app),
           str(frozen / 'src/data/h2o_ccpvdz_psi4.FCIDUMP'), str(args.iterations), str(args.interval),
           str(args.reference_energy), str(target), 'olsen-only', str(args.z_threshold),
           str(args.max_wavefunction_size), '64', '1', '0', str(args.z_threshold), '0']
    save(output / 'command.json', {'argv': cmd})
    for sig in (signal.SIGUSR1, signal.SIGTERM, signal.SIGINT):
        signal.signal(sig, stop_signal)
    deadline = time.monotonic() + args.budget_seconds - 600
    process = None
    interrupted = False
    try:
        print('START one joint solve, fixed blocks=', args.iterations, 'output=', output, flush=True)
        process = subprocess.Popen(cmd, start_new_session=True)
        last_plotted = 0
        while process.poll() is None:
            if STOP or time.monotonic() >= deadline:
                interrupted = True
                os.killpg(process.pid, signal.SIGTERM)
                try:
                    process.wait(timeout=30)
                except subprocess.TimeoutExpired:
                    os.killpg(process.pid, signal.SIGKILL)
                break
            time.sleep(5)
            state = analyze(output)
            if state['reports'] >= last_plotted + 10:
                analyze(output, plots=True)
                last_plotted = state['reports']
        process.wait()
        state = analyze(output, plots=True)
        if interrupted or process.returncode != 0 or not state['complete']:
            raise RuntimeError(f'Joint solve stopped ({process.returncode}); partial JSONL retained')
        save(output / 'completion.json', {'status': 'complete', 'iterations': args.iterations,
                                         'finished_unix': time.time()})
    except BaseException as error:
        if process is not None and process.poll() is None:
            os.killpg(process.pid, signal.SIGTERM)
            try:
                process.wait(timeout=30)
            except subprocess.TimeoutExpired:
                os.killpg(process.pid, signal.SIGKILL)
                process.wait()
        save(output / 'completion.json', {'status': 'partial_or_failed', 'reason': str(error)})
        raise
    finally:
        analyze(output, plots=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest='command', required=True)
    p = sub.add_parser('prepare')
    for name in ('root', 'frozen', 'eigen', 'reference-file'):
        p.add_argument('--' + name, type=Path, required=True)
    p = sub.add_parser('seal'); p.add_argument('--frozen', type=Path, required=True)
    p = sub.add_parser('analyze'); p.add_argument('--output', type=Path, required=True)
    p = sub.add_parser('run')
    p.add_argument('--frozen', type=Path, required=True); p.add_argument('--output', type=Path, required=True)
    p.add_argument('--iterations', type=int, default=5000000); p.add_argument('--interval', type=int, default=10000)
    p.add_argument('--budget-seconds', type=int, default=64800)
    p.add_argument('--max-wavefunction-size', type=int, default=848256040)
    p.add_argument('--z-threshold', type=float, default=0)
    p.add_argument('--reference-energy', type=float, default=-76.2418600697)
    args = parser.parse_args()
    if args.command == 'prepare': prepare_single(args)
    elif args.command == 'seal': seal(args.frozen)
    elif args.command == 'analyze': print(json.dumps(analyze(args.output, True), indent=2))
    else:
        if args.iterations <= 0 or args.interval <= 0 or args.iterations % args.interval:
            parser.error('Positive fixed iterations divisible by report interval required')
        if args.budget_seconds <= 600 or args.max_wavefunction_size <= 0:
            parser.error('Invalid budget/capacity')
        if not math.isfinite(args.z_threshold) or args.z_threshold != 0:
            parser.error('This uncompressed diagnostic requires z_threshold=0')
        if args.reference_energy != -76.2418600697:
            parser.error('This experiment uses the certified H2O reference -76.2418600697')
        run(args)


if __name__ == '__main__':
    main()
