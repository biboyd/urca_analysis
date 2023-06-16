#!/bin/bash
#SBATCH --account=m3018
#SBATCH --nodes=1
#SBATCH --tasks-per-node=8
#SBATCH --mem=48G
#SBATCH --qos=shared
#SBATCH --time=00:15:00
#SBATCH --constraint=cpu
#SBATCH --mail-user=boyd.brendan@stonybrook.edu
#SBATCH --mail-type=FAIL

export OMP_NUM_THREADS=8
export MPICH_MAX_THREAD_SAFETY=multiple
export MY_EXEC=./main3d.gnu.x86-milan.MPI.ex

fname=$1
my_inputs=$@
srun -n 8 $MY_EXEC ~/my_urca_analysis/nu_inputs problem.plotfilename="$fname" problem.outfile="${fname%/}_nu" $my_inputs