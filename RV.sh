#!/bin/bash

#PBS -l select=1:ncpus=1:mem=200m
#PBS -l walltime=10:00:00

cd /mnt/scratch/ws/svdremov/202305240321Exp_Data/RV_Eq/Damp_0.01/RV_contin_200k

make

./dzeq

echo "SCV is gotta go, sir"



