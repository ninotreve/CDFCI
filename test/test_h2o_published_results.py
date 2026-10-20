"""Audit the committed H2O data without Slurm, FCIDUMP or a solver binary."""
import csv
import hashlib
import json
from pathlib import Path
import shutil
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / 'experiment_results/energy_correction/h2o_ccpvdz/1027249'
sys.path.insert(0, str(ROOT / 'experiments/energy_correction'))
import h2o_single_run as single


def read(name):
    return json.loads((DATA / name).read_text())


@unittest.skipUnless((DATA / 'completion.json').is_file(), 'Published data absent in source-only snapshot')
class PublishedH2O(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.manifest = read('manifest.json')
        cls.summary = read('summary.json')
        cls.raw = read('runs/joint.json')
        cls.events = single.journal(DATA / 'runs/joint.json.jsonl')
        cls.rows = [r for r in cls.events if r['record_type'] == 'report']

    def test_complete_fixed_length_single_solve(self):
        self.assertEqual(read('completion.json')['status'], 'complete')
        self.assertTrue(self.summary['complete'])
        self.assertEqual(read('plan.json')['number_of_solves'], 1)
        self.assertEqual(self.summary['last_iteration'], 5000000)
        self.assertEqual(len(self.rows), 500)
        self.assertEqual([r['iteration'] for r in self.rows], list(range(10000, 5000001, 10000)))
        self.assertEqual(sum(e['record_type'] == 'run_start' for e in self.events), 1)
        self.assertEqual(sum(e['record_type'] == 'run_end' for e in self.events), 1)
        opt = self.raw['options']
        self.assertTrue(opt['energy_correction']['enabled'])
        self.assertTrue(opt['energy_correction']['olsen_enabled'])
        self.assertEqual(opt['stopping_dx_threshold'], 0)
        self.assertEqual(opt['z_threshold'], 0)
        self.assertEqual(opt['num_coordinates'], 64)
        self.assertEqual(self.raw['trajectory_run'], self.summary['solve_summary'])

    def test_raw_streamed_and_csv_energies_agree(self):
        with (DATA / 'three_way_trajectory.csv').open() as stream:
            csv_rows = list(csv.DictReader(stream))
        self.assertEqual(len(self.raw['trajectory']), 500)
        self.assertEqual(len(csv_rows), 500)
        ref = self.manifest['reference_energy']
        for raw, event, row in zip(self.raw['trajectory'], self.rows, csv_rows):
            self.assertEqual(raw['iteration'], event['iteration'])
            self.assertEqual(int(row['iteration']), event['iteration'])
            for field in ('variational_energy', 'external_correction', 'internal_correction',
                          'corrected_energy', 'correction_seconds'):
                self.assertEqual(raw[field], event[field])
            energies = (event['variational_energy'],
                        event['variational_energy'] + event['external_correction'],
                        event['corrected_energy'])
            self.assertAlmostEqual(energies[2], energies[1] + event['internal_correction'], delta=1e-12)
            self.assertTrue(ref < energies[2] < energies[1] < energies[0])
            for name, value in zip(single.ENERGIES, energies):
                self.assertEqual(float(row[name + '_energy']), value)
                self.assertEqual(float(row[name + '_signed_error']), value - ref)

    def test_reanalysis_reproduces_committed_tables_and_summary(self):
        with tempfile.TemporaryDirectory(prefix='cdfci-h2o-published-') as directory:
            output = Path(directory)
            shutil.copy2(DATA / 'manifest.json', output / 'manifest.json')
            (output / 'runs').mkdir()
            shutil.copy2(DATA / 'runs/joint.json.jsonl', output / 'runs/joint.json.jsonl')
            result = single.analyze(output)
            self.assertEqual(result, self.summary)
            for name in ('three_way_trajectory.csv', 'target_tolerances.csv'):
                self.assertEqual((output / name).read_bytes(), (DATA / name).read_bytes())

    def test_reference_and_executed_source_hashes(self):
        reference = DATA / 'reference/h2o_ccpvdz_new.out'
        self.assertEqual(hashlib.sha256(reference.read_bytes()).hexdigest(),
                         self.manifest['reference_policy']['sha256'])
        # Documentation/CMake evolve after publication; executed numerical code does not.
        for name in ('include/solver.h', 'include/wavefunction.h', 'include/energy_correction.h',
                     'test/benchmark_energy_correction.cpp',
                     'experiments/energy_correction/h2o_single_run.py',
                     'experiments/energy_correction/h2o_experiment.py'):
            self.assertEqual(hashlib.sha256((ROOT / name).read_bytes()).hexdigest(),
                             self.manifest['source']['files']['src/' + name], name)


if __name__ == '__main__':
    unittest.main()
