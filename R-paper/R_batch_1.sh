#!/bin/bash
#
#$ -cwd
#$ -V
#$ -j y
#$ -S /bin/bash
#$ -m beas
#$ -M wcai@berkeley.edu
#$ -l h=!(compute-0-9)
#

module load gcc-7.3.0
mpirun -n 1 R --vanilla < scenario_fix_t.R > scenario_fix_t.Rout

