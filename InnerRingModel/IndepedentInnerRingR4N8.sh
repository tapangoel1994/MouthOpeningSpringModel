#! /bin/bash
#SBATCH -t 2:00:00
#SBATCH -n 16
#SBATCH -N 1
#SBATCH -p compute

mkdir -p $HOME/matlabdata

mkdir $HOME/matlabdata/$SLURM_JOBID

cd $HOME/MouthOpening/NewModel/IndependentInnerRing

module load matlab

export MATLAB_PREFDIR = $HOME/matlabdata/$SLURM_JOBID

matlab -batch "run('./IndepedentInnerRingR4N8.m');exit"