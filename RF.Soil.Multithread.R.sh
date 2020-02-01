#!/bin/sh -l
#PBS -N RF.Soil.R
#PBS -q woeste
#PBS -l nodes=1:ppn=24
#PBS -l walltime=96:00:00
#PBS -m abe

module load r

cd /scratch/brown/will1809/RF
R -f RF.Soil.Multithread.R
