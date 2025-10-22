#!/bin/bash
#SBATCH --job-name=calc_eos
#SBATCH --account=m3018
#SBATCH --nodes=1
#SBATCH --tasks-per-node=8
#SBATCH --mem=64G
#SBATCH --qos=shared
#SBATCH --time=00:30:00
#SBATCH --constraint=cpu
#SBATCH --mail-user=boyd.brendan@stonybrook.edu
#SBATCH --mail-type=FAIL
#SBATCH --output=sl-eosderiv-%j.out

export MY_EOS_EXEC=feosderiv3d.gnu.x86-milan.MPI.ex


echo Doing $@ plotfiles

for fname_raw in $@
do

    # get file name
	fname=${fname_raw%/}
	base_fname=${fname#plotfiles/}

	# run each calculation
    echo Running eos on $fname
	srun -n 8 $MY_EOS_EXEC diag.plotfile=$fname

	# move to repsective directories
	mv eosderiv.${base_fname} eos_plotfiles/
done
