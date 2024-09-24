#!/bin/bash
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --mem=64G
#SBATCH --partition=short-40core
#SBATCH --time=02:00:00
#SBATCH --mail-user=boyd.brendan@stonybrook.edu
#SBATCH --mail-type=FAIL

module clear -f
module load slurm
module load python
module load openmpi/gcc12.1/4.1.4

source activate nucastro

cd /gpfs/projects/CalderGroup/BrendanSpace/plots/hires_urca_problem/radvel_volume_render

bash generate_urca_renders.sh
