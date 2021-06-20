#!/bin/bash
#SBATCH -J bc-verify
#SBATCH -o %x.%j.out
#SBATCH -p nocona
#SBATCH -N 2
#SBATCH --ntasks-per-node=128
#SBATCH -t 24:00:00
#SBATCH -D /home/zhansen/Dist_Mem_GAPBS/Dist_Mem_GAPBS/Dist_Mem_GAPBS

module load gcc/10.1.0 openmpi/4.0.4

cd /home/zhansen/Dist_Mem_GAPBS/gapbs/src
g++ -std=c++11 -o BC bc.cc

cd /home/zhansen/Dist_Mem_GAPBS/Dist_Mem_GAPBS/Dist_Mem_GAPBS
oshc++ -std=c++11 -o BC bc.cc

cd /home/zhansen/Dist_Mem_GAPBS/Dist_Mem_GAPBS/Dist_Mem_GAPBS
srun -N 2 -n 4 ./BC -u 10 -n 5 -v
mv bc_output.txt ../../bc_output.txt
/home/zhansen/Dist_Mem_GAPBS/gapbs/src/BC -u 10 -n 5 -v
echo "############# Trial Complete! #############\n"

srun -N 2 -n 16 ./BC -u 15 -n 5 -v
mv bc_output.txt ../../bc_output.txt
/home/zhansen/Dist_Mem_GAPBS/gapbs/src/BC -u 15 -n 5 -v
echo "############# Trial Complete! #############\n"

srun -N 2 -n 8 ./BC -g 15 -n 5 -v
mv bc_output.txt ../../bc_output.txt
/home/zhansen/Dist_Mem_GAPBS/gapbs/src/BC -g 15 -n 5 -v
echo "############# Trial Complete! #############\n"

srun -N 2 -n 16 ./BC -g 15 -n 5 -v
mv bc_output.txt ../../bc_output.txt
/home/zhansen/Dist_Mem_GAPBS/gapbs/src/BC -g 15 -n 5 -v
echo "############# Trial Complete! #############\n"

srun -N 2 -n 4 ./BC -f /home/zhansen/Dist_Mem_GAPBS/gapbs/benchmark/graphs/raw/twitter.el -n 5 -v
mv bc_output.txt ../../bc_output.txt
/home/zhansen/Dist_Mem_GAPBS/gapbs/src/BC -f /home/zhansen/Dist_Mem_GAPBS/gapbs/benchmark/graphs/raw/twitter.el -n 5 -v
echo "############# Trial Complete! #############\n"

#srun -N 1 -n 32 ./BC -f /home/zhansen/Dist_Mem_GAPBS/gapbs/benchmark/graphs/raw/twitter.el -n 1 -v
#mv bc_output.txt ../../bc_output.txt
#/home/zhansen/Dist_Mem_GAPBS/gapbs/src/BC -f /home/zhansen/Dist_Mem_GAPBS/gapbs/benchmark/graphs/raw/twitter.el -n 1 -v
#echo "############# Trial Complete! #############\n"
