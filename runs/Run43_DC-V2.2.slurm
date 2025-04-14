#!/bin/bash
#
#SBATCH --job-name=Run43_DC-V2.2
#SBATCH --output=Run43_DC-V2.2.log
#
#SBATCH --ntasks=200
#SBATCH --cpus-per-task=1
#SBATCH --partition=year-long-cpu

# These variables allow system load to show a decrease as processes finish their jobs
export I_MPI_THREAD_YIELD=3
export I_MPI_THREAD_SLEEP=100

# This library could be added to the path (it may still be, but set it here)
export I_MPI_PMI_LIBRARY=/usr/lib64/libpmi2.so

cd /home/jvallino/Projects/MCP-DEB/MEP-DEB/Darwin-Chemostat/V2/

module list
srun ../Darwin_Chemostat Run43_DC-V2.2
