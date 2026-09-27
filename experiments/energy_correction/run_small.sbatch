#!/bin/bash
#SBATCH -p amd_256
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH -J cdfci_pe_small
#SBATCH -t 00:30:00
#SBATCH --output=cdfci_pe_small_%j.out
#SBATCH --error=cdfci_pe_small_%j.out

set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
exec bash "${SCRIPT_DIR}/run_all.sbatch" small
