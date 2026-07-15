#!/bin/bash
#SBATCH --job-name=G_PTA_MPI
#SBATCH --output=OUTPUT/out.txt
#SBATCH --error=OUTPUT/err.txt
#SBATCH --ntasks=11
#SBATCH --cpus-per-task=1
#SBATCH --time=01:00:00

cd $SLURM_SUBMIT_DIR/SOURCE

module load mpi/openmpi-x86_64

mpirun -np 11 ./GA_par.exe