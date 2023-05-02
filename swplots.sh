#!/bin/bash

files=$@
for file in $files
do
	echo "generating plots for $file"
	# generate the plots
	python ~/my_urca_analysis/plot_figures.py $file "radial_velocity"  --zlo=-1e1 --zup=1e1 --linthresh=1e-2 --width=1.5e3 &
	python ~/my_urca_analysis/plot_figures.py $file "Hnuc" --zlo=-1e8 --zup=1e8 --linthresh=1e4 --width=1.5e3 &
    python ~/my_urca_analysis/plot_figures.py $file "Ye_asymmetry" --zlo=-2.6e-5 --zup=-1e-5 --width=1.5e3 &
    python ~/my_urca_analysis/plot_figures.py $file "X(na23)" --width=1.5e3 &
    python ~/my_urca_analysis/plot_figures.py $file "X(ne23)" --width=1.5e3 &
    python ~/my_urca_analysis/plot_figures.py $file magvel --zlo='min' --zup='max' --width=1.5e3 &

	# mv the plt files to plotfiles directory so we dont pick them up again
	wait
done

wait
