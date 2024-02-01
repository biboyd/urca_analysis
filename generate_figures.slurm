#!/bin/bash
#SBATCH --job-name=plot-urca
#SBATCH --account=m3018
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --cpus-per-task=1
#SBATCH --qos=shared
#SBATCH --time=00:25:00
#SBATCH --constraint=cpu
#SBATCH --mail-type=FAIL

export CONV_EXE=main3d.gnu.x86-milan.MPI.ex
module load cpu

mkdir -p plotfiles
mkdir -p profiles
mkdir -p conv_projs
mkdir -p trash
#mv plt?????00 trash/
rm -r plt?????00
rm -r plt?????00.old.*
wait

for file in plt???????
do
	# check if any plots to make
    if [ $file = "plt???????" ];
	then
		echo "no plots to generate"
	else
		echo "generating plots for $file"
	    # generate the plots
	    for fld in "Ye_asymmetry" "radial_velocity" "Hnuc" "MachNumber" "c12_complement"
	    do
            if [ $fld = 'c12_complement' ]
	    	then
	    		python ~/my_urca_analysis/plot_figures.py $file "$fld" -s -lc white &
	    	elif [ $fld = 'radial_velociy' ]
	    	then
	    		python ~/my_urca_analysis/plot_figures.py $file "$fld" -cH &
	    	else
	    		python ~/my_urca_analysis/plot_figures.py $file "$fld" &
	    	fi
	    done
	    
	    # generate regular projections
        python ~/my_urca_analysis/saveprofile.py $file 
		wait
		python ~/my_urca_analysis/plot_profiles.py profiles/${file}.csv $file

	    # generate convective criteria projections
        #SBATCH --mem=32G
		#wait
	    #srun -n 8 $CONV_EXE ~/my_urca_analysis/conv_inputs problem.plotfilename=$file problem.proj_outfile=conv_projs/${file}.txt
	    #python ~/my_urca_analysis/plot_conv_proj.py $file conv_projs/${file}.txt

	    # mv the plt files to plotfiles directory so we dont pick them up again
	    wait
	    echo "moving $file to plotfiles/$file"
	    mv $file "plotfiles/"
	fi
done

wait
