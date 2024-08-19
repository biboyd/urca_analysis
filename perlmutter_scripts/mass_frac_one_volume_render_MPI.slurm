#!/bin/bash
#SBATCH --account=m3018
#SBATCH --job-name=frac-render
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --mem=64G
#SBATCH --qos=shared
#SBATCH --time=00:10:00
#SBATCH --constraint=cpu
#SBATCH --mail-user=boyd.brendan@stonybrook.edu
#SBATCH --mail-type=FAIL

export OMP_NUM_THREADS=8
export MPICH_MAX_THREAD_SAFETY=multiple

my_inputs=$@
echo $my_inputs
#srun -n 4 python ~/my_urca_analysis/mass_frac_one_volume_render.py $my_inputs
srun python ~/my_urca_analysis/mass_frac_one_volume_render.py $my_inputs
