#!/usr/bin/env python3
"""Frozen, calibrated same-allocation H2O experiment; also analyze interrupted JSONL."""
from __future__ import annotations

import argparse
import csv
import hashlib
import io
import importlib
import json
import math
import os
from pathlib import Path
import shutil
import signal
import statistics
import subprocess
import time

STOP = False


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def save(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + '.tmp')
    temporary.write_text(json.dumps(value, indent=2, allow_nan=False) + '\n')
    temporary.replace(path)


def read(path):
    return json.loads(Path(path).read_text())


def csv_file(path, rows):
    if not rows:
        return
    stream = io.StringIO(newline='')
    writer = csv.DictWriter(stream, fieldnames=list(rows[0]), lineterminator='\n')
    writer.writeheader()
    writer.writerows(rows)
    temporary = path.with_name(path.name + '.tmp')
    temporary.write_text(stream.getvalue())
    temporary.replace(path)


def prepare(root, frozen, eigen):
    root, frozen = root.resolve(), frozen.resolve()
    if (frozen / 'source_manifest.json').exists():
        raise ValueError('Refusing to overwrite an existing frozen snapshot')
    names = subprocess.check_output(
        ['git', 'ls-files', '--cached', '--others', '--exclude-standard', '-z'], cwd=root).decode().split('\0')
    for name in sorted(set(names)):
        if name == 'CMakeLists.txt' or name.startswith(('include/', 'src/', 'test/', 'experiments/energy_correction/')):
            source, target = root / name, frozen / 'src' / name
            if source.is_file():
                target.parent.mkdir(parents=True, exist_ok=True)
                shutil.copy2(source, target)
    for name in ('h2o_ccpvdz_psi4.FCIDUMP', 'h2o_sto3g_psi4.FCIDUMP'):
        target = frozen / 'src/data' / name
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(root / 'data' / name, target)
    target = frozen / 'src/regression_tests/h2o/ccpvdz_psi4/output'
    target.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(root / 'regression_tests/h2o/ccpvdz_psi4/output', target)
    shutil.copytree(eigen, frozen / 'eigen', ignore=shutil.ignore_patterns('.git', '__pycache__'))
    files = {str(p.relative_to(frozen)): digest(p) for p in (frozen / 'src').rglob('*') if p.is_file()}
    save(frozen / 'source_manifest.json', {
        'git_commit': subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=root, text=True).strip(),
        'working_tree_status': subprocess.check_output(['git', 'status', '--short'], cwd=root, text=True),
        'files': files, 'eigen_core_sha256': digest(frozen / 'eigen/Eigen/Core'),
        'note': 'Copied actual working-tree inputs before submission, including uncommitted fixes; all phases use this snapshot.'})


def seal(frozen):
    manifest = read(frozen / 'source_manifest.json')
    manifest.update(binary_sha256=digest(frozen / 'build/bin/cdfci_energy_correction_benchmark_omp'),
                    compiler=subprocess.check_output(['/usr/bin/c++', '--version'], text=True).splitlines()[0],
                    build_type='Release', build_flags='-O3 -DNDEBUG',
                    sealed_at_unix=time.time())
    save(frozen / 'source_manifest.json', manifest)


def verify(frozen, all_sources=False):
    manifest = read(frozen / 'source_manifest.json')
    if digest(frozen / 'build/bin/cdfci_energy_correction_benchmark_omp') != manifest['binary_sha256']:
        raise RuntimeError('Frozen binary changed')
    if all_sources:
        for name, value in manifest['files'].items():
            if digest(frozen / name) != value:
                raise RuntimeError('Frozen input changed: ' + name)
    return manifest


def journal(path):
    """Ignore only an unterminated last line; reject corruption elsewhere."""
    if not path.exists():
        return []
    lines = path.read_text().splitlines(keepends=True)
    events = [json.loads(line) for line in lines if line.endswith('\n')]
    rows = [e for e in events if e.get('record_type') == 'report']
    if len({r['iteration'] for r in rows}) != len(rows):
        raise ValueError('Multiple/duplicate runs in a single-mode journal')
    return events


def reports(path):
    return [e for e in journal(path) if e.get('record_type') == 'report']


def plan_steps(pilots, cap, interval, pairs, available, reserve, safety=1.5):
    rates = {}
    for name, data in pilots.items():
        rows, total = data
        if len(rows) < 2:
            raise ValueError('Calibration needs at least two report points per mode')
        samples = [(b['wall_seconds'] - a['wall_seconds']) / (b['iteration'] - a['iteration'])
                   for a, b in zip(rows, rows[1:])]
        rates[name] = max(total / rows[-1]['iteration'], max(samples))
    rate = max(rates.values())
    count = 1 + 2 * pairs
    steps = min(cap, int((available - reserve) / (safety * count * rate))) // interval * interval
    if steps < 2 * interval:
        raise RuntimeError('Insufficient safe allocation budget for repeated equal-length runs')
    return {'iterations': steps, 'pairs': pairs, 'number_of_solves': count,
            'seconds_per_step_by_mode': rates, 'slowest_rate': rate, 'safety_factor': safety,
            'available_seconds_after_calibration': available, 'reserve_seconds': reserve,
            'estimated_production_seconds_with_safety': steps * rate * count * safety,
            'requested_iteration_cap': cap,
            'note': 'Worst calibration rate including allocation and interval growth, 1.5x safety and a separate time reserve; not a completion guarantee.'}


def schedule(pairs):
    entries = [('raw', 'raw')]
    for i in range(pairs):
        order = ('ip', 'joint') if i % 2 == 0 else ('joint', 'ip')
        entries.extend((f'{mode}_{i}', mode) for mode in order)
    return entries


def analyze(output, plots=False):
    config = read(output / 'manifest.json')
    selected = read(output / 'plan.json') if (output / 'plan.json').exists() else None
    entries = schedule(config['pairs'])
    available = {}
    targets = []
    ref = config['reference_energy']
    for label, mode in entries:
        events = journal(output / 'runs' / (label + '.json.jsonl'))
        rows = [r for r in events if r.get('record_type') == 'report']
        if not rows:
            continue
        start = next(e for e in events if e.get('record_type') == 'run_start')
        flags = start['options']['energy_correction']
        if flags['enabled'] != (mode != 'raw') or (mode != 'raw' and flags['olsen_enabled'] != (mode == 'joint')):
            raise ValueError('Invalid correction switches in ' + label)
        if mode == 'ip' and any(r['internal_correction'] != 0 for r in rows):
            raise ValueError('Unexpected Olsen contribution in IP-only')
        available[label] = rows
        for tolerance in (1e-3, 1e-4, 1e-5, 1e-6):
            hit = next((r for r in rows if abs(r['variational_energy' if mode == 'raw' else 'corrected_energy'] - ref) <= tolerance), None)
            targets.append({'run': label, 'mode': mode, 'tolerance_hartree': tolerance,
                            'iteration': hit['iteration'] if hit else None,
                            'wall_seconds': hit['wall_seconds'] if hit else None,
                            'hamiltonian_columns': hit['hamiltonian_columns'] if hit else None})
    save(output / 'available_trajectories.json', available)
    csv_file(output / 'target_tolerances.csv', targets)
    summaries = []
    for i in range(config['pairs']):
        ip, joint = available.get(f'ip_{i}'), available.get(f'joint_{i}')
        if not ip or not joint:
            continue
        a, b = {r['iteration']: r for r in ip}, {r['iteration']: r for r in joint}
        common = sorted(a.keys() & b.keys())
        if not common:
            continue
        step = common[-1]
        for k in common:
            if a[k]['hamiltonian_columns'] != b[k]['hamiltonian_columns']:
                raise ValueError('Hamiltonian-work mismatch')
        completed = bool(selected and step == selected['iterations'])
        ends = []
        for label in (f'ip_{i}', f'joint_{i}'):
            ev = journal(output / 'runs' / (label + '.json.jsonl'))
            ends.append(next((e['summary']['seconds'] for e in ev if e.get('record_type') == 'run_end'), None))
        use_solve = completed and all(x is not None for x in ends)
        ti, tj = ends if use_solve else (a[step]['wall_seconds'], b[step]['wall_seconds'])
        summaries.append({'pair': i, 'iterations': step, 'complete': bool(use_solve),
                          'timing_basis': 'full_solve' if use_solve else 'common_prefix_report_clock',
                          'ip_seconds': ti, 'joint_seconds': tj,
                          'olsen_overhead_vs_ip': tj / ti - 1})
    complete = [s for s in summaries if s['complete']]
    values = [s['olsen_overhead_vs_ip'] for s in complete]
    summary = {'planned_iterations': selected['iterations'] if selected else None,
               'observed_reports_by_run': {label: len(rows) for label, rows in available.items()},
               'planned_pairs': config['pairs'], 'complete_pairs': len(complete),
               'paired_timings': summaries,
               'olsen_overhead_median': statistics.median(values) if values else None,
               'olsen_overhead_min': min(values) if values else None,
               'olsen_overhead_max': max(values) if values else None,
               'note': 'Paired medians/range are descriptive, not confidence intervals; raw measured once. Same allocation, balanced alternating pair order. Partial pairs are not pooled into full-solve statistics.',
               'reference_policy': config['reference_policy']}
    save(output / 'summary.json', summary)
    if not all(k in available for k in ('raw', 'ip_0', 'joint_0')):
        return summary
    maps = {k: {r['iteration']: r for r in available[k]} for k in ('raw', 'ip_0', 'joint_0')}
    common = sorted(set(maps['raw']) & set(maps['ip_0']) & set(maps['joint_0']))
    rows = []
    for k in common:
        row = {'iteration': k}
        for label, mode in [('raw', 'raw'), ('ip_0', 'ip'), ('joint_0', 'joint')]:
            r = maps[label][k]
            row.update({mode + '_energy': r['variational_energy' if mode == 'raw' else 'corrected_energy'],
                        mode + '_wall_seconds': r['wall_seconds'],
                        mode + '_stored_determinants': r['stored_determinants'],
                        mode + '_stored_wavefunction_entries': r['stored_wavefunction_entries'],
                        mode + '_hamiltonian_columns': r['hamiltonian_columns']})
            row[mode + '_signed_error'] = row[mode + '_energy'] - ref
        row['joint_report_seconds'] = maps['joint_0'][k]['correction_seconds']
        rows.append(row)
    csv_file(output / 'three_way_trajectory.csv', rows)
    if not plots or not rows:
        return summary
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    floor = config['reference_policy']['plot_resolution_guard_hartree']
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.4))
    for mode, label in [('raw', 'raw CDFCI'), ('ip', 'IP only (pair 0)'), ('joint', 'IP + Olsen (pair 0)')]:
        resolved = [r for r in rows if abs(r[mode + '_signed_error']) >= floor]
        for ax, key in [(axes[0], '_stored_determinants'), (axes[1], '_wall_seconds')]:
            ax.loglog([r[mode + key] for r in resolved], [abs(r[mode + '_signed_error']) for r in resolved], label=label)
            ax.set_ylim(bottom=floor)
    axes[0].set(xlabel='own nonzero variational coefficients', ylabel='absolute error (Ha)')
    axes[1].set(xlabel='own wall time at report (s)', ylabel='absolute error (Ha)')
    accumulated = []
    total = 0
    for r in rows:
        total += r['joint_report_seconds']
        accumulated.append(total / r['joint_wall_seconds'])
    axes[2].semilogx([r['joint_stored_wavefunction_entries'] for r in rows], accumulated)
    axes[2].set(xlabel='joint stored wavefunction entries', ylabel='direct report time / joint elapsed time')
    for ax in axes:
        ax.grid(True, which='both', alpha=.25)
    for ax in axes[:2]: ax.legend()
    fig.suptitle('Representative pair 0; reference-resolution guard applied; direct timer is not Olsen marginal cost')
    fig.tight_layout()
    for ext in ('pdf', 'png'): fig.savefig(output / ('three_way_advantage.' + ext), dpi=180)
    plt.close(fig)
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))
    for mode in ('raw', 'ip', 'joint'):
        axes[0].plot(common, [r[mode + '_signed_error'] for r in rows], label=mode)
    axes[0].axhspan(-floor, floor, color='gray', alpha=.2, label='unresolved diagnostic region')
    axes[0].set_yscale('symlog', linthresh=floor)
    axes[0].set(xlabel='block iteration', ylabel='signed E - E_ref (Ha)')
    axes[0].legend()
    if summaries:
        axes[1].bar([str(s['pair']) + ('' if s['complete'] else '*') for s in summaries],
                    [100 * s['olsen_overhead_vs_ip'] for s in summaries])
    axes[1].axhline(0, color='black', linewidth=.8)
    axes[1].set(xlabel='IP/Olsen pair (* = partial prefix)', ylabel='(T_joint - T_IP) / T_IP (%)')
    fig.tight_layout()
    for ext in ('pdf', 'png'): fig.savefig(output / ('signed_errors_and_olsen_overhead.' + ext), dpi=180)
    plt.close(fig)
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))
    for mode in ('raw', 'ip', 'joint'):
        selected_rows = [r for r in rows if abs(r[mode + '_signed_error']) >= floor]
        axes[0].semilogy([r['iteration'] for r in selected_rows], [abs(r[mode + '_signed_error']) for r in selected_rows], label=mode)
        axes[1].loglog([r[mode + '_hamiltonian_columns'] for r in selected_rows], [abs(r[mode + '_signed_error']) for r in selected_rows], label=mode)
    axes[0].set(xlabel='block iteration', ylabel='absolute error (Ha)', ylim=(floor, None))
    axes[1].set(xlabel='actual Hamiltonian columns', ylabel='absolute error (Ha)', ylim=(floor, None))
    for ax in axes:
        ax.legend(); ax.grid(True, which='both', alpha=.25)
    fig.tight_layout()
    for ext in ('pdf', 'png'): fig.savefig(output / ('errors_vs_iteration_and_work.' + ext), dpi=180)
    plt.close(fig)
    return summary


def stop_signal(signum, frame):
    global STOP
    STOP = True


def run(args):
    global STOP
    frozen, output = args.frozen.resolve(), args.output.resolve()
    source = verify(frozen, True)
    # Fail before calibration, rather than discovering missing plotting packages
    # after spending most of a production allocation.
    for dependency in ('numpy', 'matplotlib'):
        importlib.import_module(dependency)
    output.mkdir(parents=True, exist_ok=True)
    if (output / 'manifest.json').exists():
        raise RuntimeError('Refusing to overwrite an existing experiment')
    config = vars(args).copy()
    config = {k: str(v) if isinstance(v, Path) else v for k, v in config.items()}
    config.update(source=source, job_id=os.environ.get('SLURM_JOB_ID'), node=os.environ.get('SLURMD_NODENAME'),
                  partition=os.environ.get('SLURM_JOB_PARTITION'), omp_threads=os.environ.get('OMP_NUM_THREADS'),
                  batch_process_affinity_cpu_ids=sorted(os.sched_getaffinity(0)),
                  allocated_cpus_per_node=os.environ.get('SLURM_JOB_CPUS_PER_NODE'),
                  binding_note='64 OpenMP threads, not a claim of 64 physical cores; each srun prints its binding mask.',
                  omp_places=os.environ.get('OMP_PLACES'), omp_proc_bind=os.environ.get('OMP_PROC_BIND'),
                  reference_policy={'reference_energy': args.reference_energy,
                     'source': 'frozen regression_tests/h2o/ccpvdz_psi4/output' if args.reference_energy == -76.2418600697 else 'user override, not certified by stored regression',
                     'rounding_half_unit_hartree': 5e-11 if args.reference_energy == -76.2418600697 else None,
                     'plot_resolution_guard_hartree': args.reference_guard,
                     'note': 'Conservative display guard, not a rigorous error bound. Negative errors inside the shaded region do not establish failure; zero printed tail drift does not certify accuracy.'},
                  zero_threshold_note='z_threshold=0 is an uncompressed-residual diagnostic, not a compressed production performance claim.')
    save(output / 'manifest.json', config)
    allocation_start = time.time()
    deadline = allocation_start + args.budget_seconds - 600
    signal.signal(signal.SIGUSR1, stop_signal)
    signal.signal(signal.SIGTERM, stop_signal)
    signal.signal(signal.SIGINT, stop_signal)
    app = frozen / 'build/bin/cdfci_energy_correction_benchmark_omp'
    fcidump = frozen / 'src/data/h2o_ccpvdz_psi4.FCIDUMP'
    modes = {'raw': 'raw-only', 'ip': 'ip-only', 'joint': 'olsen-only'}

    def execute(label, mode, steps, area):
        verify(frozen)
        path = output / area / (label + '.json')
        path.parent.mkdir(parents=True, exist_ok=True)
        if path.exists() or Path(str(path) + '.jsonl').exists():
            raise RuntimeError('Refusing to overwrite a run journal: ' + str(path))
        cmd = ['srun', '--ntasks=1', '--cpus-per-task=64', '--cpu-bind=verbose,cores', str(app), str(fcidump),
               str(steps), str(min(args.interval, steps)), str(args.reference_energy), str(path), modes[mode],
               str(args.z_threshold), str(args.max_wavefunction_size), str(args.coordinates), '1', '0', str(args.z_threshold), '0']
        print('START', area, label, 'mode', mode, 'steps', steps, flush=True)
        start = time.monotonic()
        process = subprocess.Popen(cmd, start_new_session=True)
        interrupted = False
        while process.poll() is None:
            if STOP or time.time() >= deadline:
                interrupted = True
                os.killpg(process.pid, signal.SIGTERM)
                try: process.wait(timeout=30)
                except subprocess.TimeoutExpired:
                    os.killpg(process.pid, signal.SIGKILL)
                break
            time.sleep(5)
            # Cheap, full-precision snapshots survive before the final JSON.
            analyze(output, plots=False)
        process.wait()
        analyze(output, plots=True)
        if interrupted or process.returncode != 0:
            raise RuntimeError(f'{label} stopped ({process.returncode}); journals and common-prefix analysis retained')
        rows = reports(Path(str(path) + '.jsonl'))
        if not rows or rows[-1]['iteration'] != steps:
            raise RuntimeError('Run did not reach the requested fixed length')
        return rows, time.monotonic() - start

    try:
        pilots = {}
        for mode in ('raw', 'ip', 'joint'):
            execute('calibration_' + mode, 'raw', 200, 'warmups')
            pilots[mode] = execute(mode, mode, args.pilot_steps, 'calibration')
        selected = plan_steps(pilots, args.max_iterations, args.interval, args.pairs,
                              args.budget_seconds - (time.time() - allocation_start), args.reserve_seconds)
        save(output / 'plan.json', selected)
        print('CALIBRATED PLAN', json.dumps(selected), flush=True)
        for label, mode in schedule(args.pairs):
            execute('production_' + label, 'raw', 200, 'warmups')
            execute(label, mode, selected['iterations'], 'runs')
        save(output / 'completion.json', {'status': 'complete', 'finished_unix': time.time()})
    except Exception as error:
        save(output / 'completion.json', {'status': 'partial_or_failed', 'reason': str(error)})
        raise
    finally:
        analyze(output, plots=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest='command', required=True)
    p = sub.add_parser('prepare'); p.add_argument('--root', type=Path, required=True); p.add_argument('--frozen', type=Path, required=True); p.add_argument('--eigen', type=Path, required=True)
    p = sub.add_parser('seal'); p.add_argument('--frozen', type=Path, required=True)
    p = sub.add_parser('analyze'); p.add_argument('--output', type=Path, required=True)
    p = sub.add_parser('run'); p.add_argument('--frozen', type=Path, required=True); p.add_argument('--output', type=Path, required=True)
    p.add_argument('--max-iterations', type=int, default=4000000); p.add_argument('--interval', type=int, default=10000)
    p.add_argument('--pilot-steps', type=int, default=50000); p.add_argument('--pairs', type=int, default=4)
    p.add_argument('--budget-seconds', type=int, default=64800); p.add_argument('--reserve-seconds', type=int, default=7200)
    p.add_argument('--coordinates', type=int, default=64); p.add_argument('--max-wavefunction-size', type=int, default=848256040)
    p.add_argument('--z-threshold', type=float, default=0); p.add_argument('--reference-energy', type=float, default=-76.2418600697)
    p.add_argument('--reference-guard', type=float, default=1e-7)
    args = parser.parse_args()
    if args.command == 'prepare': prepare(args.root, args.frozen, args.eigen)
    elif args.command == 'seal': seal(args.frozen)
    elif args.command == 'analyze': print(json.dumps(analyze(args.output, True), indent=2))
    else:
        if args.pairs < 2 or args.pairs % 2 or args.interval <= 0 or args.pilot_steps < 2 * args.interval or args.pilot_steps % args.interval:
            raise ValueError('Use an even number of >=2 pairs and >=2 full calibration reports')
        if args.max_iterations < 2 * args.interval or args.max_iterations % args.interval or args.coordinates != 64:
            raise ValueError('Invalid iteration cap or coordinate count')
        if not all(math.isfinite(x) for x in (args.reference_energy, args.reference_guard, args.z_threshold)) or args.reference_guard <= 0 or args.z_threshold < 0:
            raise ValueError('Invalid finite reference/threshold settings')
        run(args)


if __name__ == '__main__':
    main()
