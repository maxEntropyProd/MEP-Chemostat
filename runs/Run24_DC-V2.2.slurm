#!/bin/bash
#
#SBATCH --job-name=Run24_DC-V2.2
#SBATCH --output=Run24_DC-V2.2.log
#
#SBATCH --ntasks=200
#SBATCH --cpus-per-task=1
#SBATCH --partition=year-long-cpu

# These variables allow system load to show a decrease as processes finish their jobs
export I_MPI_THREAD_YIELD=3
export I_MPI_THREAD_SLEEP=100

# This library could be added to the path (it may still be, but set it here)
export I_MPI_PMI_LIBRARY=/usr/lib64/libpmi2.so

module list
srun ../src/Darwin_Chemostat Run24_DC-V2.2
