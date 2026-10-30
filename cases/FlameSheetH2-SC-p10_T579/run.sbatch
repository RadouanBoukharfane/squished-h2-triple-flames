#!/bin/bash
#SBATCH --job-name=pele_1d_SC_p10_T579
#SBATCH --partition=<your_partition>
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --time=03:00:00
#SBATCH --account=<your_account>
#SBATCH --qos=<your_qos>
#SBATCH --output=slurm-%x-%j.out
# 1D flame phi = 1.6, p10_T579, constant-Schmidt transport (constant Schmidt numbers)
module purge
module load <your MPI/compiler environment>  # e.g. foss/2024a
ulimit -s unlimited
cd "$SLURM_SUBMIT_DIR"
hostname; date
mpirun -np $SLURM_NTASKS ./PeleLMeX2d.gnu.MPI.SanDiegoH2.ex input.1d_phi1.6 > pele.log 2>&1
echo "exit=$?"; date
