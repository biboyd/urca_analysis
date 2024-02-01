#!/bin/bash
#SBATCH --job-name=profiles
#SBATCH --account=m3018
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --qos=shared
#SBATCH --time=00:45:00
#SBATCH --constraint=cpu

mkdir -p profiles

files=$@

for f in $files
do
	python ~/my_urca_analysis/save_all_profile.py $f
done
