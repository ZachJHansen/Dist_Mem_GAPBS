#!/bin/bash
#SBATCH -J tc-verify
#SBATCH -o %x.%j.out
#SBATCH -p nocona
#SBATCH -N 2
#SBATCH --ntasks-per-node=128
#SBATCH -t 36:00:00
#SBATCH -D /home/zhansen/Dist_Mem_GAPBS/Dist_Mem_GAPBS/Dist_Mem_GAPBS

module load gcc/10.1.0 openmpi/4.0.4

cd /home/zhansen/Dist_Mem_GAPBS/gapbs/src
g++ -std=c++11 -o TC tc.cc

cd /home/zhansen/Dist_Mem_GAPBS/Dist_Mem_GAPBS/Dist_Mem_GAPBS
oshc++ -std=c++11 -o TC tc.cc

cd /home/zhansen/Dist_Mem_GAPBS/Dist_Mem_GAPBS/Dist_Mem_GAPBS
srun -N 2 -n 4 ./TC -u 10 -n 3 -v
mv tc_output.txt ../../tc_output.txt
/home/zhansen/Dist_Mem_GAPBS/gapbs/src/TC -u 10 -n 3 -v
echo "############# Trial Complete! #############\n"

srun -N 2 -n 16 ./TC -u 15 -n 3 -v
mv tc_output.txt ../../tc_output.txt
/home/zhansen/Dist_Mem_GAPBS/gapbs/src/TC -u 15 -n 3 -v
echo "############# Trial Complete! #############\n"

srun -N 2 -n 8 ./TC -g 15 -n 3 -v
mv tc_output.txt ../../tc_output.txt
/home/zhansen/Dist_Mem_GAPBS/gapbs/src/TC -g 15 -n 3 -v
echo "############# Trial Complete! #############\n"

srun -N 2 -n 16 ./TC -g 15 -n 3 -v
mv tc_output.txt ../../tc_output.txt
/home/zhansen/Dist_Mem_GAPBS/gapbs/src/TC -g 15 -n 3 -v
echo "############# Trial Complete! #############\n"

srun -N 2 -n 4 ./TC -sf /home/zhansen/Dist_Mem_GAPBS/gapbs/benchmark/graphs/raw/twitter.el -n 5 -v
mv tc_output.txt ../../tc_output.txt
/home/zhansen/Dist_Mem_GAPBS/gapbs/src/TC -sf /home/zhansen/Dist_Mem_GAPBS/gapbs/benchmark/graphs/raw/twitter.el -n 5 -v
echo "############# Trial Complete! #############\n"

#srun -N 1 -n 32 ./TC -f /home/zhansen/Dist_Mem_GAPBS/gapbs/benchmark/graphs/raw/twitter.el -n 1 -v
#mv tc_output.txt ../../tc_output.txt
#/home/zhansen/Dist_Mem_GAPBS/gapbs/src/TC -f /home/zhansen/Dist_Mem_GAPBS/gapbs/benchmark/graphs/raw/twitter.el -n 1 -v
#echo "############# Trial Complete! #############\n"
