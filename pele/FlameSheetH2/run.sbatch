#!/bin/bash
#SBATCH --job-name=pele_1d
#SBATCH --partition=<your_partition>
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --time=03:00:00
#SBATCH --account=<your_account>
#SBATCH --qos=<your_qos>
#SBATCH --output=slurm-%x-%j.out
# PeleLMeX verification: 1D H2/air flame phi = 1.6, 1 atm, 300 K, San Diego + mixture-averaged (vs Cantera).
module purge
module load <your MPI/compiler environment>  # e.g. foss/2024a
ulimit -s unlimited
cd "$SLURM_SUBMIT_DIR"
hostname; date
mpirun -np $SLURM_NTASKS ./PeleLMeX2d.gnu.MPI.SanDiegoH2.ex input.1d_phi1.6 > pele.log 2>&1
echo "exit=$?"; date
