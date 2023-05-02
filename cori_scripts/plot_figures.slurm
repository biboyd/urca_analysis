#!/bin/bash
#SBATCH --job-name=plot-urca
#SBATCH --account=m3018
#SBATCH --nodes=1
#SBATCH --ntasks=3
#SBATCH --cpus-per-task=2
#SBATCH --qos=shared
#SBATCH --time=00:05:00
#SBATCH --constraint=haswell

mkdir -p plotfiles

for file in plt???????
do
	echo "generating plots for $file"
	# generate the plots
	for fld in "Ye_asymmetry" "radial_velocity" "Hnuc" "MachNumber" "c12_complement" "mu" "A23_frac"
	do
		python ~/my_urca_analysis/plot_figures.py $file "$fld" &
	done

	# mv the plt files to plotfiles directory so we dont pick them up again
	wait
	echo "moving $file to plotfiles/$file"
	mv $file "plotfiles/"
done

wait
