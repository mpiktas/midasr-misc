#!/bin/sh
#SBATCH -p long
#SBATCH -n21 
R_LIBS=/scratch/lustre/zemlys/lib64/R/library
export R_LIBS

mpirun  -n 1 /soft/zemlys/R/bin/R --vanilla --slave -f /scratch/lustre/zemlys/R/midasr-misc/jvcheck/20clsim.R
