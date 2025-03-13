#!/bin/sh
##PJM settings
#PJM -L "rscunit=ito-a"
#PJM -L "rscgrp=ito-ss"
#PJM -L "vnode=1"
#PJM -L "vnode-core=36"
#PJM -L "elapse=96:00:00"
#GE settings
#$ -N job_name
#$ -pe smp 32
#$ -cwd
#$ -V
#$ -q all.q
#$ -S /bin/bash
#Torque settings
##PBS -j oe
##PBS -o out.o
#PBS -l nodes=1:ppn=16
#Slurm settings
#SBATCH -p salmon
#SBATCH -n 1
#SBATCH -c 16
#SBATCH -J job_name
#SBATCH -o out.o
#SBATCH -e out.o

#export OMP_NUM_THREADS=32
#cd $PBS_O_WORKDIR
python main.py
