#!/bin/bash
#
#$ -cwd
#$ -V
#$ -j y
#$ -S /bin/bash
#

# module load gcc-4.9.4
module load gcc-7.3.0
# mpirun -n 1 R --vanilla < scenario.R > scenario.Rout
R CMD BATCH --vanilla < scenario.R > scenario.Rout

