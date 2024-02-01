#!/bin/bash
#SBATCH --job-name=plot-sw-urca
#SBATCH --account=m3018
#SBATCH --nodes=1
#SBATCH --ntasks=6
#SBATCH --cpus-per-task=1
#SBATCH --qos=shared
#SBATCH --time=00:15:00
#SBATCH --constraint=cpu

module load cpu

mkdir -p plotfiles
rm -r plt?????00

for file in plt???????
do

	# check if any plots to make
    if [ $file = "plt???????" ];
	then
		echo "no plots to generate"
    else
        echo "generating plots for $file"
        # generate the plots
        bash ~/my_urca_analysis/swplots.sh $file
        # mv the plt files to plotfiles directory so we dont pick them up again
        wait
        echo "moving $file to plotfiles/$file"
        mv $file "plotfiles/"
    fi
done

wait
