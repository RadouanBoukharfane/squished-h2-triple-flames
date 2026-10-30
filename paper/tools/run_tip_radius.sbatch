#!/bin/bash
#SBATCH --job-name=tip_radius
#SBATCH --partition=<your_partition>
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --time=02:00:00
#SBATCH --account=<your_account>
#SBATCH --qos=<your_qos>
#SBATCH --output=slurm-%x-%j.out
# Figure 4 post-processing: tip radius of the PeleLMeX triple flames (tools/tip_radius.py).
cd "$SLURM_SUBMIT_DIR"; date
for TT in 1000 1500 2000; do T_TIP=$TT python -u tools/tip_radius.py data/tip_radius_T$TT.json \
  C1 C4 C5 C6 A6P-erf A6P-tanh B1-p10-D20 B1-p10-D5b B1-p10-D5 B1-p10-D3b B1-p20-D20 B1-p20-D7 B1-p20-D5 B1-p20-D3; echo "== T=$TT done"; done
echo "exit=$?"; date
