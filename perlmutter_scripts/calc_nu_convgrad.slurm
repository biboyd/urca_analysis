#!/bin/bash
#SBATCH --job-name=nuloss-convgrad
#SBATCH --account=m3018
#SBATCH --nodes=1
#SBATCH --tasks-per-node=8
#SBATCH --mem=48G
#SBATCH --qos=shared
#SBATCH --time=00:30:00
#SBATCH --constraint=cpu
#SBATCH --mail-user=boyd.brendan@stonybrook.edu
#SBATCH --mail-type=FAIL
#SBATCH --output=sl-nuloss-convgrad-%j.out

export MY_NU_EXEC=fneutrinos3d.gnu.x86-milan.MPI.ex
export MY_CONV_EXEC=fconvgrad3d.gnu.x86-milan.MPI.ex


echo Doing $@ plotfiles

for fname_raw in $@
do

    # get file name
	fname=${fname_raw%/}
	base_fname=${fname#plotfiles/}

	# run each calculation
    echo Running Nuloss on $fname
	srun -n 8 $MY_NU_EXEC diag.plotfile=$fname
    echo Running convgrad on $fname
	srun -n 8 $MY_CONV_EXEC diag.plotfile=$fname diag.spherical=1

	# move to repsective directories
    echo Moving nu_loss.${base_fname}  and convgrad.${base_fname} 
	mv nu_loss.${base_fname} nuloss_plotfiles/
	mv convgrad.${base_fname} convgrad_plotfiles/
done
