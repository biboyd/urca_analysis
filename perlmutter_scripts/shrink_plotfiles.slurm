#!/bin/bash
#SBATCH --job-name=shrink_pfile
#SBATCH --account=m3018
#SBATCH --nodes=1
#SBATCH --tasks-per-node=8
#SBATCH --mem=48G
#SBATCH --qos=shared
#SBATCH --time=00:30:00
#SBATCH --constraint=cpu
#SBATCH --mail-user=boyd.brendan@stonybrook.edu
#SBATCH --mail-type=FAIL
#SBATCH --output=sl-shrink-%j.out

export MY_EXEC=fshrink3d.gnu.x86-milan.MPI.ex


echo Doing $@ plotfiles

for fname_raw in $@
do
    # get file name
	fname=${fname_raw%/}
	base_fname=${fname#plotfiles/}

	# run each calculation
    echo Shrinking $fname
	srun -n 8 $MY_EXEC diag.plotfile=$fname fab.format=NATIVE_32

	# move to directory and copy basestate/header info
	mv small.${base_fname} small_plotfiles/
    cp ${fname_raw}/job_info ${fname_raw}/Base?C_0 small_plotfiles/small.${base_fname}/
done
