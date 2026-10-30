#!/bin/bash
#SBATCH --job-name=cantera_sweep
#SBATCH --partition=<your_partition>
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=56
#SBATCH --time=04:00:00
#SBATCH --account=<your_account>
#SBATCH --qos=<your_qos>
#SBATCH --output=slurm-%x-%j.out
# Pressure sweep of the consumption Markstein number along the isentropic path (markstein_sweep.py).
cd "$SLURM_SUBMIT_DIR"
hostname; date
python -u markstein_sweep.py markstein_sweep.json 9
echo "exit=$?"; date
