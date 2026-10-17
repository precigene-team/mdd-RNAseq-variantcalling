#!/bin/bash
# [[04_asm/07_drach_analysis.sbatch]]
# Goal: Run DRACH motif discovery in consensus peaks and intersect with per-sample filtered SNPs.

#SBATCH --job-name=drach_analysis
#SBATCH --output=_data/logs/drach/%x_%j.log
#SBATCH --error=_data/logs/drach/%x_%j.err
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=02:00:00

set -e
set -o pipefail

if [ -n "${SLURM_SUBMIT_DIR:-}" ]; then
    export PROJECT_ROOT="$SLURM_SUBMIT_DIR"
else
    export PROJECT_ROOT="${PROJECT_ROOT:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
fi

eval "$(micromamba shell hook --shell bash 2>/dev/null || mamba shell hook --shell bash 2>/dev/null || conda shell.bash hook 2>/dev/null)"
micromamba activate mdd 2>/dev/null || mamba activate mdd 2>/dev/null || conda activate mdd

mkdir -p "$PROJECT_ROOT/_data/logs/drach"

echo ">>> Step 04: Finding DRACH motif sites in consensus peaks..."
python3 "$PROJECT_ROOT/04_asm/04_find_drach_sites.py"

echo ">>> Step 06: Intersecting filtered SNPs with DRACH sites and annotating disruption..."
bash "$PROJECT_ROOT/04_asm/06_intersect_drach_snps.sh"

echo ">>> DRACH motif analysis complete!"
