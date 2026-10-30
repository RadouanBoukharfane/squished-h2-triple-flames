#!/bin/bash
#SBATCH --job-name=pele_COLDF-B1-p20-D4-SC
#SBATCH --partition=<your_partition>
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=56
#SBATCH --time=06:00:00
#SBATCH --exclusive
#SBATCH --account=<your_account>
#SBATCH --qos=<your_qos>
#SBATCH --output=slurm-%x-%j.out
# PeleLMeX reference run of case C1 (Chaplet et al. 2026).
module purge
module load <your MPI/compiler environment>  # e.g. foss/2024a
ulimit -s unlimited
cd "$SLURM_SUBMIT_DIR"
hostname; date
mpirun -np $SLURM_NTASKS ./PeleLMeX2d.gnu.MPI.SanDiegoH2.ex input.case > pele.log 2>&1
echo "exit=$?"; date
