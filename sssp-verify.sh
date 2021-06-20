#!/bin/bash
#SBATCH -J sssp-verify
#SBATCH -o %x.%j.out
#SBATCH -p nocona
#SBATCH -N 2
#SBATCH --ntasks-per-node=128
#SBATCH -t 24:00:00
#SBATCH -D /home/zhansen/Dist_Mem_GAPBS/Dist_Mem_GAPBS/Dist_Mem_GAPBS

module load gcc/10.1.0 openmpi/4.0.4

cd /home/zhansen/Dist_Mem_GAPBS/gapbs/src
g++ -std=c++11 -o SSSP sssp.cc

cd /home/zhansen/Dist_Mem_GAPBS/Dist_Mem_GAPBS/Dist_Mem_GAPBS
oshc++ -std=c++11 -o SSSP new_sssp.cc

cd /home/zhansen/Dist_Mem_GAPBS/Dist_Mem_GAPBS/Dist_Mem_GAPBS
srun -N 2 -n 4 ./SSSP -u 10 -n 3 -v
mv sssp_output.txt ../../sssp_output.txt
/home/zhansen/Dist_Mem_GAPBS/gapbs/src/SSSP -u 10 -n 3 -v
echo "############# Trial Complete! #############\n"

srun -N 2 -n 16 ./SSSP -u 15 -n 3 -v
mv sssp_output.txt ../../sssp_output.txt
/home/zhansen/Dist_Mem_GAPBS/gapbs/src/SSSP -u 15 -n 3 -v
echo "############# Trial Complete! #############\n"

srun -N 2 -n 8 ./SSSP -g 15 -n 3 -v
mv sssp_output.txt ../../sssp_output.txt
/home/zhansen/Dist_Mem_GAPBS/gapbs/src/SSSP -g 15 -n 3 -v
echo "############# Trial Complete! #############\n"

srun -N 2 -n 16 ./SSSP -g 15 -n 3 -v
mv sssp_output.txt ../../sssp_output.txt
/home/zhansen/Dist_Mem_GAPBS/gapbs/src/SSSP -g 15 -n 3 -v
echo "############# Trial Complete! #############\n"

srun -N 2 -n 4 ./SSSP -f /home/zhansen/Dist_Mem_GAPBS/gapbs/benchmark/graphs/raw/twitter.el -n 5 -v
mv sssp_output.txt ../../sssp_output.txt
/home/zhansen/Dist_Mem_GAPBS/gapbs/src/SSSP -f /home/zhansen/Dist_Mem_GAPBS/gapbs/benchmark/graphs/raw/twitter.el -n 5 -v
echo "############# Trial Complete! #############\n"

#srun -N 1 -n 32 ./TC -f /home/zhansen/Dist_Mem_GAPBS/gapbs/benchmark/graphs/raw/twitter.el -n 1 -v
#mv sssp_output.txt ../../sssp_output.txt
#/home/zhansen/Dist_Mem_GAPBS/gapbs/src/TC -f /home/zhansen/Dist_Mem_GAPBS/gapbs/benchmark/graphs/raw/twitter.el -n 1 -v
#echo "############# Trial Complete! #############\n"
