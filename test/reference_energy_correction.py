"""Optional independent PySCF FCI reference for a benchmark FCIDUMP."""
import json
import os
import sys
from pathlib import Path

import pyscf
from pyscf import fci, lib
from pyscf.tools import fcidump

lib.num_threads(int(os.environ.get("PYSCF_NUM_THREADS", "2")))
source, destination = map(Path, sys.argv[1:3])
# PySCF's comma-joined parser rejects a line containing only '&FCI'.
# Join that line to the following header line; preserve every integral.
compatible = destination.with_suffix(".FCIDUMP")
compatible.write_text(source.read_text().replace("&FCI\n", "&FCI ", 1))
data = fcidump.read(str(compatible))
solver = fci.direct_spin1.FCI()
solver.conv_tol = 1e-13
solver.max_cycle = 200
energy, vector = solver.kernel(
    data["H1"], data["H2"], data["NORB"], data["NELEC"], ecore=data["ECORE"]
)
result = dict(energy=energy, converged=bool(solver.converged), norb=data["NORB"],
              nelec=data["NELEC"], conv_tol=solver.conv_tol, pyscf=pyscf.__version__,
              source=str(source))
destination.write_text(json.dumps(result, indent=2) + "\n")
print(json.dumps(result))
if not solver.converged:
    raise SystemExit("FCI reference did not converge")
