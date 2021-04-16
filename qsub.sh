#!/bin/bash
#SBATCH -J qsub_tc
#SBATCH -o %x.%j.out
#SBATCH -p nocona
#SBATCH -N 2
#SBATCH --ntasks-per-node=128
#SBATCH -t 24:00:00
#SBATCH -D /home/zhansen/Dist_Mem_GAPBS

module load gcc/10.1.0 openmpi/4.0.4
oshc++ -std=c++11 -o TC tc.cc

oshrun --npernode 32 ./TC -g 25 -n 1 -v
echo "Trial 1 Complete!"
