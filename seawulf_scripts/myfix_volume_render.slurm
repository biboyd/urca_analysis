#!/bin/bash
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --mem=64G
#SBATCH --partition=short-28core
#SBATCH --time=00:15:00
#SBATCH --mail-user=boyd.brendan@stonybrook.edu
#SBATCH --mail-type=FAIL


source activate my-env

my_inputs=$@
echo $my_inputs
srun -n 1 python ~/Repo/urca_analysis/volume_rendering/myfix_volume_render.py $my_inputs
