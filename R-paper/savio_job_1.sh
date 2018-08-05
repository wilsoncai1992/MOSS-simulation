#!/bin/bash
#SBATCH --job-name=MOSS
#SBATCH -A co_biostat
#SBATCH -p savio2
#SBATCH -n 72
#SBATCH --exclusive
#SBATCH -t 24:00:00
#SBATCH --mail-user=wcai@berkeley.edu
#SBATCH --output=slurm.out
#SBATCH --error=slurm.out

# module load gcc
module load r/3.4.2
module load r-packages
module load ml/superlearner/current-r-3.4.2
module load rmpi/0.6-6
mpirun R CMD BATCH --no-save scenario.R scenario.Rout