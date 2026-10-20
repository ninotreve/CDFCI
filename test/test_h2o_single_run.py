"""Single-run scheduling, signed curves and recovery, without Slurm access."""
import argparse
import csv
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'experiments/energy_correction'))
import h2o_single_run as single


def events(steps, reference=-2., complete=False):
    result = [{'record_type': 'run_start', 'options': {
        'num_iterations': steps, 'stopping_dx_threshold': 0, 'z_threshold': 0,
        'energy_correction': {'enabled': True, 'olsen_enabled': True}}}]
    for step in range(10000, steps + 1, 10000):
        result.append({'record_type': 'report', 'iteration': step, 'phase': 'ip+olsen',
                       'wall_seconds': step / 1000, 'variational_energy': reference + .5,
                       'external_correction': -.3, 'internal_correction': -.25,
                       'corrected_energy': reference - .05, 'stored_determinants': step,
                       'stored_wavefunction_entries': step * 2,
                       'hamiltonian_columns': step * 64 - 63, 'correction_seconds': .01})
    if complete:
        result.append({'record_type': 'run_end', 'summary': {'iterations': steps, 'seconds': 123.}})
    return result


def write_events(path, data):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(''.join(json.dumps(e) + '\n' for e in data))


class SingleRun(unittest.TestCase):
    def test_same_state_signed_errors_and_partial_recovery(self):
        with tempfile.TemporaryDirectory() as directory:
            out = Path(directory)
            single.save(out / 'manifest.json', {'iterations': 5000000, 'z_threshold': 0,
                        'reference_energy': -2., 'reference_policy': {'source': 'test'}})
            path = out / 'runs/joint.json.jsonl'
            data = events(20000)
            data[0]['options']['num_iterations'] = 5000000
            write_events(path, data)
            with path.open('a') as stream:
                stream.write('{"unfinished":')
            state = single.analyze(out)
            self.assertFalse(state['complete'])
            self.assertEqual(state['last_iteration'], 20000)
            self.assertEqual(state['below_reference_reports']['ip_olsen'], 2)
            self.assertEqual(state['below_reference_reports']['ip'], 0)
            with (out / 'three_way_trajectory.csv').open() as stream:
                row = next(csv.DictReader(stream))
            self.assertAlmostEqual(float(row['variational_signed_error']), .5)
            self.assertAlmostEqual(float(row['ip_signed_error']), .2)
            self.assertAlmostEqual(float(row['ip_olsen_signed_error']), -.05)
            self.assertEqual(float(row['wall_seconds_joint']), 10.)
            self.assertNotIn('raw_wall_seconds', row)

    def test_exactly_one_five_million_step_joint_launch(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            frozen, out = root / 'frozen', root / 'output'
            reference = frozen / 'src/reference/h2o_ccpvdz_new.out'
            reference.parent.mkdir(parents=True)
            reference.write_text('certified reference fixture\n')
            args = argparse.Namespace(command='run', frozen=frozen, output=out, iterations=5000000,
                                      interval=10000, budget_seconds=64800,
                                      max_wavefunction_size=848256040, z_threshold=0,
                                      reference_energy=-76.2418600697)
            commands = []

            class FinishedProcess:
                returncode = 0
                def __init__(self, command, **kwargs):
                    commands.append(command)
                    path = Path(command[9] + '.jsonl')
                    write_events(path, events(int(command[6]), float(command[8]), complete=True))
                def poll(self): return 0
                def wait(self, **kwargs): return 0

            with patch.object(single, 'verify', return_value={}), \
                    patch.object(single.subprocess, 'Popen', FinishedProcess), \
                    patch.object(single, 'make_plots'), patch.object(single.signal, 'signal'):
                single.run(args)
            self.assertEqual(len(commands), 1)
            self.assertEqual(commands[0][6], '5000000')
            self.assertEqual(commands[0][10], 'olsen-only')
            self.assertEqual(commands[0][-1], '0')
            self.assertEqual(single.read(out / 'plan.json')['number_of_solves'], 1)
            self.assertFalse(single.read(out / 'plan.json')['automatic_iteration_reduction'])
            self.assertEqual(single.read(out / 'completion.json')['status'], 'complete')
            self.assertEqual(single.read(out / 'summary.json')['reports'], 500)
            self.assertFalse((out / 'calibration').exists())
            self.assertFalse((out / 'warmups').exists())
            self.assertEqual(len(list((out / 'runs').glob('*.jsonl'))), 1)


if __name__ == '__main__':
    unittest.main()
