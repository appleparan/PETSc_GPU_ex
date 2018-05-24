#!/bin/bash

#SBATCH -J PETSc_TEST                   # Job name
#SBATCH -o jobout.%j   # Name of stdout output file (%j expands to %jobId)
#SBATCH -n 1
#SBATCH -p gpu                           # queue or partiton name
#SBATCH --gres=gpu:1                      # Num Devices

mpiexec -N 1 $PWD/poisson
# End of File.
