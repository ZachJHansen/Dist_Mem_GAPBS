#!/bin/bash
#SBATCH -J TC_Synthetic
#SBATCH -o %x.o%j
#SBATCH -p nocona
#SBATCH -N 1 
#SBATCH --ntasks-per-node=128
#SBATCH -t 4:00:00
#SBATCH -D /home/zhansen/Dist_Mem_GAPBS
#SBATCH --mail-user=zhansen@ttu.edu
#SBATCH --mail-type=ALL
 
module load gcc/10.1.0 openmpi/4.0.4
oshc++ -std=c++11 -o TC tc.cc
oshrun -N 6 ./TC -g 20 -n 8

