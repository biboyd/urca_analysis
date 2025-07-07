#!/bin/bash

for p in plt0*/
do 
	python ~/my_urca_analysis/save_plot_profiles/save_all_profile.py $p
    for axis in x z
    do   
        # radial velocity
    	python ~/my_urca_analysis/plot_figures.py $p  radial_velocity -a $axis -w 2e3 --zlo='-30' --zup='30' --linthresh=0.5 -o plots_radvel_wide/ -sph 170 412 490
    
        # vorticity
        python ~/my_urca_analysis/plot_figures.py $p  vort -a $axis -w 2e3 --zlo=0. -g False -sph 170 412 490
    
        # Nuc Energy
        python ~/my_urca_analysis/plot_figures.py $p  Hnuc  -a $axis -w 2e3 --zlo='-1e12' --zup='1e12' --linthresh=1e8 -sph 170 412 490
    
        # composition
    	for fld in X\(Mg24\) A23_ratio X\(Ne23\) X\(Na23\) X\(C12\)
    	do 
    		python ~/my_urca_analysis/plot_figures.py $p "$fld" -g False -a $axis  -w 2e3 -o select_comp_slices/ -sph 170 412 490
    	done
    done
    
    mv $p plotfiles/
done

python ~/my_urca_analysis/overtime_calculations/convect_overtime.py
