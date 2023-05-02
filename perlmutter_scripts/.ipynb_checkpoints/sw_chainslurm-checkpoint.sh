#!/bin/sh -f


if [ ! "$1" ]; then
  echo "usage: chainslurm.sh jobid number script"
  echo "       set jobid -1 for no initial dependency"
  exit -1
fi

if [ ! "$2" ]; then
  echo "usage: chainslurm.sh jobid number script"
  echo "       set jobid -1 for no initial dependency"
  exit -1
fi

if [ ! "$3" ]; then
  echo "usage: chainslurm.sh jobid number script"
  echo "       set jobid -1 for no initial dependency"
  exit -1
fi

module load cpu

oldjob=$1
numjobs=$2
script=$3

if [ $numjobs -gt "20" ]; then
    echo "too many jobs requested"
    exit -1
fi

firstcount=1

if [ $oldjob -eq "-1" ]
then 
    echo chaining $numjobs jobs
    
    echo starting job 1 with no dependency
    aout=`sbatch --parsable ${script}`
    echo "   " jobid: $aout
    echo " "
    sleep 1
    echo "and adding a plot script onto job 1"
    sbatch -d afterany:${aout} ~/my_urca_analysis/plot_sw_figures.slurm
	echo " "
    oldjob=$aout
    firstcount=2
    sleep 3
else
    echo chaining $numjobs jobs starting with $oldjob
fi

for count in `seq $firstcount 1 $numjobs`
do
  echo starting job $count to depend on $oldjob
  aout=`sbatch --parsable -d afterany:${oldjob} ${script}`
  echo "   " jobid: $aout
  echo " "
  sleep 1
  echo "and adding a plot script onto job $count"
  sbatch -d afterany:${aout} ~/my_urca_analysis/plot_sw_figures.slurm
  echo " "
  oldjob=$aout
  sleep 1
done
