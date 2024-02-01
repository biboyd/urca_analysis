#!/bin/bash
#SBATCH --account=m3018
#SBATCH --nodes=1
#SBATCH --tasks-per-node=16
#SBATCH --mem=32G
#SBATCH --qos=shared
#SBATCH --time=00:15:00
#SBATCH --constraint=cpu
#SBATCH --mail-user=boyd.brendan@stonybrook.edu
#SBATCH --mail-type=FAIL

export OMP_NUM_THREADS=4
export MPICH_MAX_THREAD_SAFETY=multiple

my_inputs=$@
srun -n 1 python ~/Repo/MAESTROeX/Exec/science/urca/analysis/scripts/volume-plot-rad_vel.py $my_inputs
