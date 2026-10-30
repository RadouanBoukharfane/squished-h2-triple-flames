#!/bin/bash
#SBATCH --job-name=eps_measured
#SBATCH --partition=<your_partition>
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --time=02:00:00
#SBATCH --account=<your_account>
#SBATCH --qos=<your_qos>
#SBATCH --output=slurm-%x-%j.out
# eps* measured on the cold runs (tools/eps_measured.py). Args: OUT.json COLD-RUN ...
module purge; module load <your MPI/compiler environment>  # e.g. foss/2024a
cd "$SLURM_SUBMIT_DIR"; date
python -u tools/eps_measured.py "$@"
echo "exit=$?"; date
