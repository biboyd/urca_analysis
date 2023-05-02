#!/bin/bash
module load esslurm
nfiles=$((`ls -d chk???????/ | wc -l` -2))

allfiles=`ls -d chk???????/`

tlo=`echo ${allfiles:3:7} | sed 's/^0*//'`

diff=$((100 * $nfiles))
tup=$((($tlo) + $diff))
bkp_files=`seq -f "chk%07g" $tlo 100 $tup`

outfile=chk${tlo}_to_${tup}.tar

echo "backing up "
echo $bkp_files
echo "to $outfile"

sbatch -J "chk2048" ~/my_urca_analysis/archive_chk.slurm my_urca2048/$outfile "$bkp_files"
module unload esslurm
