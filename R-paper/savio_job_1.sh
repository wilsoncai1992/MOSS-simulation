#!/bin/bash
#SBATCH --job-name=moss
#SBATCH -A co_biostat
##SBATCH -A fc_tmlesimu
#SBATCH -p savio2
##SBATCH -p savio2_bigmem
##SBATCH --qos=savio_lowprio
#SBATCH --nodes 18
#SBATCH --mem-per-cpu 6G
#SBATCH --exclusive
#SBATCH -t 12:00:00
#SBATCH --mail-type=all
#SBATCH --mail-user=wcai@berkeley.edu
#SBATCH --output=slurm.out
#SBATCH --error=slurm.out
#SBATCH --exclude=n0088.savio2,n0089.savio2,n0098.savio2,n0131.savio2,n0292.savio2,n0190.savio2,n0291.savio2,n0092.savio2,n0104.savio2,n0048.savio2,n0144.savio2,n0190.savio2,n0135.savio2,n0139.savio2,n0041.savio2,n0044.savio2,n0046.savio2,n0076.savio2,n0201.savio2,n0033.savio2,n0083.savio2,n0064.savio2,n0055.savio2,n0035.savio2,n0027.savio2,n0028.savio2,n0037.savio2,n0038.savio2,n0039.savio2,n0064.savio2,n0059.savio2,n0036.savio2,n0051.savio2,n0053.savio2,n0061.savio2,n0058.savio2,n0050.savio2,n0124.savio2,n0136.savio2

module load r/3.5.1
module load r-packages
module load openmpi
mpirun R CMD BATCH --no-save scenario_fix_t.R
