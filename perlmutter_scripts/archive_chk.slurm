#!/bin/bash
#SBATCH -q xfer
#SBATCH -t 04:00:00

outfile=$1

chkfiles=$2

htar -cf $outfile $chkfiles

mv $chkfiles checkfiles/
