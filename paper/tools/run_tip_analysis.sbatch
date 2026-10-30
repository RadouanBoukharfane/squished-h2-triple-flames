#!/bin/bash
#SBATCH --job-name=tip_analysis
#SBATCH --partition=<your_partition>
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=56
#SBATCH --time=03:00:00
#SBATCH --account=<your_account>
#SBATCH --qos=<your_qos>
#SBATCH --output=slurm-%x-%j.out
# Review points 1 and 5: tip Markstein number and eps* seen by the tip (tools/tip_analysis.py). Args: OUT NPROC CASES...
cd "$SLURM_SUBMIT_DIR"; date
module purge; module load <your MPI/compiler environment>  # e.g. foss/2024a
python -u tools/tip_analysis.py "$@"
echo "exit=$?"; date
