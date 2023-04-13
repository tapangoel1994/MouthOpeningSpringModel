#!/bin/bash
#SBATCH -t 8:00:00
#SBATCH -n 16
#SBATCH -N 1
#SBATCH -p unowned
#SBATCH --mail-type=BEGIN,END,FAIL,REQUEUE
#SBATCH --mail-user=tgoel36@gatech.edu

mkdir -p $HOME/matlabdata

mkdir $HOME/matlabdata/$SLURM_JOBID

module load matlab

export MATLAB_PREFDIR = $HOME/matlabdata/$SLURM_JOBID

matlab -batch "run('./IndepedentOuterRingR4N8Chimera2_5NF.m');exit"