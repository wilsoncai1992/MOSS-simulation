#!/bin/bash
#SBATCH --job-name=test
#SBATCH -A co_biostat
#SBATCH -p savio2
#SBATCH -n 72
#SBATCH --exclusive
#SBATCH -t 24:00:00
#SBATCH --mail-user=wcai@berkeley.edu
#SBATCH --output=slurm.out
#SBATCH --error=slurm.out

module load r/3.2.5
mpirun R CMD BATCH --no-save scenario6.R scenario6.Rout