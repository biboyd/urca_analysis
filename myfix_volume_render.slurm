#!/bin/bash
#SBATCH --account=m3018
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --mem=64G
#SBATCH --qos=shared
#SBATCH --time=00:15:00
#SBATCH --constraint=cpu
#SBATCH --mail-user=boyd.brendan@stonybrook.edu
#SBATCH --mail-type=BEGIN,END,FAIL

export OMP_NUM_THREADS=4
export MPICH_MAX_THREAD_SAFETY=multiple

my_inputs=$@
echo $my_inputs
srun -n 1 python ~/my_urca_analysis/myfix_volume_render.py $my_inputs
