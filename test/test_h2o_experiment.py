"""Small, deterministic tests for calibration and interrupted-run analysis."""
import importlib.util
import csv
import json
from pathlib import Path
import tempfile
import unittest

spec = importlib.util.spec_from_file_location('h2o', Path(__file__).resolve().parents[1] /
    'experiments/energy_correction/h2o_experiment.py')
h2o = importlib.util.module_from_spec(spec)
spec.loader.exec_module(h2o)


class Experiment(unittest.TestCase):
    def test_balanced_order(self):
        pairs = [x[1] for x in h2o.schedule(4)][1:]
        self.assertEqual(pairs, ['ip', 'joint', 'joint', 'ip', 'ip', 'joint', 'joint', 'ip'])

    def test_budget_counts_all_nine_solves(self):
        rows = [{'iteration': 10000, 'wall_seconds': 50}, {'iteration': 20000, 'wall_seconds': 110}]
        p = h2o.plan_steps({m: (rows, 120) for m in ['raw', 'ip', 'joint']},
                          4000000, 10000, 4, 63000, 7200)
        self.assertEqual(p['number_of_solves'], 9)
        self.assertLess(p['iterations'], 4000000)
        self.assertLessEqual(p['estimated_production_seconds_with_safety'] + 7200, 63000)

    def test_incomplete_final_line_only(self):
        with tempfile.TemporaryDirectory() as d:
            path = Path(d) / 'x.jsonl'
            path.write_text('{"record_type":"metadata"}\n{"incomplete":')
            self.assertEqual(len(h2o.journal(path)), 1)
            path.write_text('{broken}\n{"incomplete":')
            with self.assertRaises(json.JSONDecodeError): h2o.journal(path)

    def test_independent_raw_and_partial_pairs(self):
        with tempfile.TemporaryDirectory() as d:
            out = Path(d)
            h2o.save(out / 'manifest.json', {'pairs': 4, 'reference_energy': -2.,
                    'reference_policy': {'plot_resolution_guard_hartree': 1e-7}})
            h2o.save(out / 'plan.json', {'iterations': 20000})
            (out / 'runs').mkdir()
            for label, mode, energy, seconds in [('raw', 'raw', -1., 10), ('ip_0', 'ip', -1.5, 20), ('joint_0', 'joint', -1.75, 21)]:
                flags = {'enabled': mode != 'raw', 'olsen_enabled': mode == 'joint'}
                row = {'record_type': 'report', 'iteration': 10000, 'wall_seconds': seconds,
                       'variational_energy': energy, 'corrected_energy': energy,
                       'stored_determinants': seconds, 'stored_wavefunction_entries': 100,
                       'hamiltonian_columns': 10000, 'internal_correction': 0,
                       'correction_seconds': .1}
                events = [{'record_type': 'run_start', 'options': {'energy_correction': flags}}, row]
                (out / 'runs' / (label + '.json.jsonl')).write_text(''.join(json.dumps(e) + '\n' for e in events))
            s = h2o.analyze(out)
            self.assertEqual(s['complete_pairs'], 0)
            self.assertAlmostEqual(s['paired_timings'][0]['olsen_overhead_vs_ip'], .05)
            self.assertIsNone(s['olsen_overhead_median'])
            text = (out / 'three_way_trajectory.csv').read_text()
            self.assertIn('raw_energy', text)
            self.assertIn('-1.0', text)
            with (out / 'three_way_trajectory.csv').open() as stream:
                combined = next(csv.DictReader(stream))
            self.assertEqual(float(combined['raw_energy']), -1.)
            self.assertEqual(float(combined['raw_wall_seconds']), 10.)
            self.assertEqual(int(combined['raw_stored_determinants']), 10)
            # A completed pair is counted only when BOTH run-end records exist.
            h2o.save(out / 'plan.json', {'iterations': 10000})
            s = h2o.analyze(out)
            self.assertEqual(s['complete_pairs'], 0)
            for label, seconds in [('ip_0', 22.), ('joint_0', 23.)]:
                with (out / 'runs' / (label + '.json.jsonl')).open('a') as stream:
                    stream.write(json.dumps({'record_type': 'run_end', 'summary': {'seconds': seconds}}) + '\n')
            s = h2o.analyze(out)
            self.assertEqual(s['complete_pairs'], 1)
            self.assertAlmostEqual(s['olsen_overhead_median'], 23 / 22 - 1)
            self.assertEqual(s['paired_timings'][0]['timing_basis'], 'full_solve')


if __name__ == '__main__': unittest.main()
