#!/bin/bash

for p in plt0*/
do 
	python ~/my_urca_analysis/save_plot_profiles/save_all_profile.py $p
	python ~/my_urca_analysis/plot_figures.py $p  radial_velocity -w 2e3 --zlo='-30' --zup='30' --linthresh=0.5 -o plots_radvel_wide/
	for fld in X\(Mg24\) A23_ratio 
	do 
		python ~/my_urca_analysis/plot_figures.py $p "$fld" -g False -a z  -w 2e3 -o select_comp_slices/
	done
	python ~/my_urca_analysis/plot_figures.py $p  Hnuc  -w 2e3 --zlo='-1e12' --zup='1e12' --linthresh=1e8 -a z
done

mv plt0*/ plotfiles/

python ~/my_urca_analysis/overtime_calculations/convect_overtime.py
