"""Check H2O script root discovery without compiling or submitting jobs."""
import os
from pathlib import Path
import shutil
import subprocess
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / 'experiments/energy_correction/run_h2o_ccpvdz.sbatch'


class BatchPaths(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(prefix='cdfci-h2o-batch-paths-')
        self.addCleanup(self.temp.cleanup)
        self.spool = Path(self.temp.name) / 'slurm' / 'job123'
        self.spool.mkdir(parents=True)
        self.copy = self.spool / 'slurm_script'
        shutil.copy2(SCRIPT, self.copy)
        self.environment = {k: v for k, v in os.environ.items()
                            if not k.startswith(('CDFCI_', 'SLURM_'))}

    def run_script(self, path, **environment):
        return subprocess.run(['bash', str(path), '--check-only'],
                              env=dict(self.environment, **environment),
                              cwd=self.temp.name, capture_output=True, text=True,
                              timeout=10)

    def assert_production_paths(self, result):
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn(f'project={ROOT}', result.stdout)
        self.assertIn(f'fcidump={ROOT}/data/h2o_ccpvdz_psi4.FCIDUMP', result.stdout)
        self.assertIn('iterations=5000000 interval=10000 z_threshold=0', result.stdout)
        self.assertIn('solves=1 mode=IP+Olsen calibration=0 timing_pairs=0', result.stdout)
        self.assertIn('fixed length; no automatic reduction', result.stdout)
        self.assertIn('threads=64 coordinates=64 max_wavefunction_size=848256040', result.stdout)
        self.assertIn('no build or solver was started', result.stdout)

    def test_slurm_spool_uses_submit_directory(self):
        self.assert_production_paths(self.run_script(self.copy, SLURM_SUBMIT_DIR=str(ROOT)))

    def test_explicit_root_overrides_submit_directory(self):
        self.assert_production_paths(self.run_script(self.copy,
                                    CDFCI_PROJECT_DIR=str(ROOT),
                                    SLURM_SUBMIT_DIR=str(self.spool)))

    def test_direct_execution_discovers_repository(self):
        self.assert_production_paths(self.run_script(SCRIPT))

    def test_invalid_submit_directory_fails_early(self):
        result = self.run_script(self.copy, SLURM_SUBMIT_DIR=str(self.spool))
        self.assertEqual(result.returncode, 2)
        self.assertIn('Invalid CDFCI project directory', result.stderr)
        self.assertIn('set CDFCI_PROJECT_DIR', result.stderr)


if __name__ == '__main__':
    unittest.main()
