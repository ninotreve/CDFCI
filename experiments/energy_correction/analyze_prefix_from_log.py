#!/usr/bin/env python3
"""Analyze a complete same-allocation prefix when IP-only final JSON is missing.

Keep stdout's rounding explicit; never fabricate an IP report timer or a
full-solve summary from a partial, timed-out invocation.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def log_rows(text: str, marker: str) -> dict[int, dict]:
    if text.count(marker) != 1:
        raise ValueError(f'Expected one unambiguous benchmark section: {marker}')
    section = text.split(marker, 1)[1].split('[Benchmark]', 1)[0]
    rows = {}
    for line in section.splitlines():
        fields = line.split()
        if len(fields) != 9 or not fields[0].isdigit():
            continue
        iteration = int(fields[0])
        if iteration in rows:
            raise ValueError(f'Duplicate report {iteration}')
        values = [float(fields[i]) for i in (1, 2, 3, 4, 8)]
        if not all(math.isfinite(v) for v in values):
            raise ValueError('Nonfinite stdout record')
        rows[iteration] = dict(iteration=iteration, variational_energy=values[0],
                               external_correction=values[1], internal_correction=values[2],
                               corrected_energy=values[3], wall_seconds=values[4],
                               stored_determinants=int(fields[6]),
                               stored_wavefunction_entries=int(fields[7]))
    return rows


def write_csv(path: Path, rows: list[dict]) -> None:
    with path.open('w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('paired_json', type=Path)
    parser.add_argument('stdout_log', type=Path)
    parser.add_argument('output', type=Path)
    parser.add_argument('--iterations', type=int, default=3000000)
    parser.add_argument('--export-sources', action='store_true',
                        help='Publish only selected JSON rows and verbatim raw/IP stdout rows, then analyze those self-contained inputs.')
    args = parser.parse_args()
    paired = json.loads(args.paired_json.read_text())
    options = paired['options']
    interval, original_steps = options['report_interval'], options['num_iterations']
    if args.iterations <= 0 or args.iterations % interval or args.iterations > original_steps:
        raise ValueError('Prefix must end at an available report boundary')
    text = args.stdout_log.read_text()
    if args.export_sources:
        args.output.mkdir(parents=True, exist_ok=True)
        json_path = args.output / 'source_paired_prefix.json'
        log_path = args.output / 'source_stdout_prefix.out'
        if json_path.exists() or log_path.exists():
            raise ValueError('Refusing to overwrite exported source files')
        selected = {k: paired[k] for k in
                    ('options', 'reference_energy', 'fcidump', 'electrons', 'spin_orbitals', 'ms2')}
        selected['trajectory'] = [r for r in paired['trajectory'] if r['iteration'] <= args.iterations]
        selected['source_export'] = {
            'analysis_iterations': args.iterations,
            'original_requested_iterations': original_steps,
            'original_paired_sha256': hashlib.sha256(args.paired_json.read_bytes()).hexdigest(),
            'original_stdout_sha256': hashlib.sha256(args.stdout_log.read_bytes()).hexdigest(),
            'note': 'Selected prefix only. Original requested iterations and stdout section markers are preserved as provenance; full-solve summaries and later reports are deliberately omitted.',
        }
        excerpts = ['# Selected same-job stdout prefix; not the complete scheduler log.',
                    f'# Only report iterations <= {args.iterations}; original invocation markers preserved.']
        for mode in ('raw CDFCI', 'IP/PT2'):
            marker = f'[Benchmark] {mode} iterations={original_steps}\n'
            if text.count(marker) != 1:
                raise ValueError('Ambiguous stdout section during export')
            section = text.split(marker, 1)[1].split('[Benchmark]', 1)[0]
            excerpts.append(marker.rstrip('\n'))
            for line in section.splitlines():
                fields = line.split()
                if line.lstrip().startswith('Iteration') or (
                        len(fields) == 9 and fields[0].isdigit()
                        and int(fields[0]) <= args.iterations):
                    excerpts.append(line)
        json_path.write_text(json.dumps(selected, indent=2) + '\n')
        log_path.write_text('\n'.join(excerpts) + '\n')
        paired, text = selected, log_path.read_text()
        args.paired_json, args.stdout_log = json_path, log_path
    raw = log_rows(text, f'[Benchmark] raw CDFCI iterations={original_steps}\n')
    ip = log_rows(text, f'[Benchmark] IP/PT2 iterations={original_steps}\n')
    joint = {r['iteration']: r for r in paired['trajectory']}
    reference = paired['reference_energy']
    rows = []
    cumulative_report_seconds = 0.0
    for iteration in range(interval, args.iterations + 1, interval):
        a, b, c = raw[iteration], ip[iteration], joint[iteration]
        if b['internal_correction'] != 0 or a['external_correction'] != 0:
            raise ValueError('Unexpected correction in raw or IP-only log')
        if abs(a['variational_energy'] - c['baseline_variational_energy']) > 5.1e-11:
            raise ValueError('Raw stdout energy does not match independent baseline JSON')
        if abs(a['wall_seconds'] - c['raw_wall_seconds']) > 0.0051:
            raise ValueError('Raw stdout time does not match baseline JSON')
        if abs(b['corrected_energy'] - b['variational_energy'] - b['external_correction']) > 1.1e-10:
            raise ValueError('IP corrected energy inconsistent with printed components')
        if c['status'] not in {'ok', 'unchecked'}:
            raise ValueError('Invalid joint correction record')
        cumulative_report_seconds += c['correction_seconds']
        row = {'iteration': iteration,
               'raw_stored_determinants': a['stored_determinants'],
               'ip_stored_determinants': b['stored_determinants'],
               'joint_stored_determinants': c['stored_determinants'],
               'joint_stored_wavefunction_entries': c['stored_wavefunction_entries'],
               'raw_energy': c['baseline_variational_energy'],
               'ip_variational_energy': b['variational_energy'],
               'ip_energy': b['corrected_energy'], 'joint_energy': c['corrected_energy'],
               'raw_wall_seconds': c['raw_wall_seconds'],
               'ip_wall_seconds': b['wall_seconds'],
               'joint_wall_seconds': c['corrected_wall_seconds'],
               'joint_cumulative_report_seconds': cumulative_report_seconds,
               'joint_report_time_fraction': cumulative_report_seconds / c['corrected_wall_seconds'],
               'olsen_overhead_vs_ip': c['corrected_wall_seconds'] / b['wall_seconds'] - 1}
        for variant in ('raw', 'ip', 'joint'):
            row[f'{variant}_error'] = abs(row[f'{variant}_energy'] - reference)
        rows.append(row)
    for key in ('raw_wall_seconds', 'ip_wall_seconds', 'joint_wall_seconds'):
        if not all(b[key] > a[key] for a, b in zip(rows, rows[1:])):
            raise ValueError('Nonmonotonic report clock')
    final = rows[-1]
    targets = []
    for tolerance in (1e-4, 1e-5, 1e-6):
        hits = {v: next((r for r in rows if r[f'{v}_error'] <= tolerance), None)
                for v in ('raw', 'ip', 'joint')}
        target = {'tolerance_hartree': tolerance}
        for variant, hit in hits.items():
            target[f'{variant}_iteration'] = hit['iteration'] if hit else None
            target[f'{variant}_seconds'] = hit[f'{variant}_wall_seconds'] if hit else None
        for variant in ('raw', 'ip'):
            target[f'{variant}_to_joint_speedup'] = (
                target[f'{variant}_seconds'] / target['joint_seconds']
                if hits[variant] and hits['joint'] else None)
        targets.append(target)
    summary = {
        'analysis_mode': 'same-allocation-complete-prefix-from-json-and-stdout',
        'source_job': '1022419', 'node': 'bigMem5',
        'original_job_state': 'TIMEOUT', 'original_requested_iterations': original_steps,
        'analysis_iterations': args.iterations, 'trajectory_points': len(rows),
        'reference_energy': reference, 'plot_resolution_floor_hartree': 1e-7,
        'wall_seconds_at_prefix_report': {v: final[f'{v}_wall_seconds'] for v in ('raw', 'ip', 'joint')},
        'olsen_overhead_vs_ip': final['olsen_overhead_vs_ip'],
        'joint_report_seconds': cumulative_report_seconds,
        'ip_report_seconds': None,
        'timing_note': 'All three trajectories are from job 1022419 on bigMem5; use their same-iteration report clocks, not full-process solve timers. One measurement per mode; fixed order and shared-node load remain confounders.',
        'log_precision': {'energy_decimal_places': 10, 'wall_seconds_decimal_places': 2,
                          'energy_rounding_half_unit_hartree': 5e-11,
                          'wall_time_rounding_half_unit_seconds': 0.005},
        'data_note': 'IP-only recovered from flushed stdout, not a completed solver JSON. Missing IP report timers and refresh counters remain unknown. Raw energies/times and joint energies/times come from original full-precision JSON; raw determinant counts come from stdout.',
        'target_tolerances': targets,
        'source_sha256': {str(p.resolve()): hashlib.sha256(p.read_bytes()).hexdigest()
                          for p in (args.paired_json, args.stdout_log)},
        'source_export': paired.get('source_export'),
    }
    args.output.mkdir(parents=True, exist_ok=True)
    (args.output / 'prefix_summary.json').write_text(json.dumps(summary, indent=2) + '\n')
    (args.output / 'prefix_trajectory.json').write_text(json.dumps(rows, indent=2) + '\n')
    write_csv(args.output / 'prefix_trajectory.csv', rows)
    write_csv(args.output / 'target_tolerances.csv', targets)
    fig, axes = plt.subplots(1, 3, figsize=(15.6, 4.4))
    for variant, label in [('raw', 'raw CDFCI'), ('ip', 'IP only'), ('joint', 'IP + Olsen')]:
        resolved = [r for r in rows if r[f'{variant}_error'] >= 1e-7]
        axes[0].loglog([r[f'{variant}_stored_determinants'] for r in resolved],
                       [r[f'{variant}_error'] for r in resolved], label=label)
        axes[1].loglog([r[f'{variant}_wall_seconds'] for r in resolved],
                       [r[f'{variant}_error'] for r in resolved], label=label)
    axes[0].set(xlabel='nonzero variational coefficients', ylabel='absolute energy error (Ha)',
                title='A. Error vs determinants', ylim=(1e-7, None))
    axes[1].set(xlabel='wall time at report (s)', ylabel='absolute energy error (Ha)',
                title='B. Error vs wall time', ylim=(1e-7, None))
    for ax in axes[:2]:
        ax.legend()
    axes[2].semilogx([r['joint_stored_wavefunction_entries'] for r in rows],
                    [r['joint_report_time_fraction'] for r in rows])
    axes[2].set(xlabel='stored wavefunction entries', ylabel='cumulative correction report time / joint time',
                title='C. Direct report-time fraction')
    for ax in axes:
        ax.grid(True, which='both', alpha=0.25)
    fig.suptitle(f'C2/cc-pVDZ: same bigMem5 allocation; prefix {args.iterations:,} steps; error >= 1e-7 Ha')
    fig.tight_layout()
    for extension in ('pdf', 'png'):
        fig.savefig(args.output / f'05_three_way_reference_resolved.{extension}', dpi=220)
    plt.close(fig)
    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.2))
    for variant, label in [('ip', 'IP only'), ('joint', 'IP + Olsen')]:
        axes[0].plot([r['iteration'] for r in rows], [r[f'{variant}_wall_seconds'] for r in rows], label=label)
    axes[0].set(xlabel='iteration', ylabel='cumulative wall time at report (s)')
    axes[0].legend()
    axes[1].semilogx([r['joint_stored_determinants'] for r in rows],
                    [100 * r['olsen_overhead_vs_ip'] for r in rows])
    axes[1].set(xlabel='nonzero variational coefficients', ylabel='extra time relative to IP (%)')
    for ax in axes:
        ax.grid(True, alpha=0.25)
    fig.suptitle('C2: same-allocation IP/Olsen comparison; IP times recovered from stdout')
    fig.tight_layout()
    for extension in ('pdf', 'png'):
        fig.savefig(args.output / f'06_ip_olsen_timing.{extension}', dpi=220)
    plt.close(fig)
    print(json.dumps(summary, indent=2))


if __name__ == '__main__':
    main()
