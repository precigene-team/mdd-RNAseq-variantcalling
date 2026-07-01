#!/bin/bash
# [[_hpc/hpc_run.sbatch]]
# Goal: Slurm submission script to execute all pipeline modules via runner.
#SBATCH --job-name=MDD_RNAseq_Prod_Run
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=48:00:00
#SBATCH --output=_data/logs/hpc_run/hpc_run_log_%j.txt
#SBATCH --error=_data/logs/hpc_run/hpc_run_error_%j.txt

# 1. Initialize Conda/Mamba (Safe hook initialization)
eval "$(micromamba shell hook --shell bash 2>/dev/null || mamba shell hook --shell bash 2>/dev/null || conda shell.bash hook 2>/dev/null)"

# 2. Activate the production environment
micromamba activate mdd 2>/dev/null || mamba activate mdd 2>/dev/null || conda activate mdd

# 3. RUN THE PIPELINE (Everything: Alignment, DGE, etc.)
./run core

echo "========================================================="
echo "   MDD RNA-seq Production Run COMPLETED SUCCESSFULLY"
echo "========================================================="
