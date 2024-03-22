#!/bin/bash
#SBATCH --nodes=1
#SBATCH --tasks-per-node=40
#SBATCH --partition=short-40core
#SBATCH --time=00:15:00
#SBATCH --mail-user=boyd.brendan@stonybrook.edu
#SBATCH --mail-type=FAIL

module load slurm
module load python
module load openmpi/gcc12.1/4.1.4
module load gnu-parallel
source activate $HOME/python-venv/nucastro/bin/activate

which python

my_wd=`pwd`
files=$@
echo $files

parallel --wd $my_wd -j 10 mpirun -np 4 python ~/Repo/urca_analysis/dipole_radvel_slice.py {} ::: $files 
