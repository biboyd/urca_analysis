#!/bin/bash
#SBATCH --account=m3018
#SBATCH --nodes=1
#SBATCH --tasks-per-node=8
#SBATCH --mem=32G
#SBATCH --qos=shared
#SBATCH --time=00:15:00
#SBATCH --constraint=cpu
#SBATCH --mail-user=boyd.brendan@stonybrook.edu
#SBATCH --mail-type=FAIL

export OMP_NUM_THREADS=8
export MPICH_MAX_THREAD_SAFETY=multiple
export MY_EXEC=./main3d.gnu.x86-milan.MPI.ex

my_inputs=$@
srun -n 8 $MY_EXEC ~/my_urca_analysis/conv_inputs $my_inputs
