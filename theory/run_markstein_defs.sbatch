#!/bin/bash
#SBATCH --job-name=cantera_markstein
#SBATCH --partition=<your_partition>
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=56
#SBATCH --time=08:00:00
#SBATCH --account=<your_account>
#SBATCH --qos=<your_qos>
#SBATCH --output=slurm-%x-%j.out
# Mechanism sensitivity and Markstein definitions (markstein_defs.py), 10 Cantera tasks in parallel.
cd "$SLURM_SUBMIT_DIR"
hostname; date
python -u markstein_defs.py markstein_defs.json 10
echo "exit=$?"; date
