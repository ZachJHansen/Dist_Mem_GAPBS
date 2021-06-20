#!/bin/bash
#SBATCH -J cc-verify
#SBATCH -o %x.%j.out
#SBATCH -p nocona
#SBATCH -N 2
#SBATCH --ntasks-per-node=128
#SBATCH -t 2:00:00
#SBATCH -D /home/zhansen/Dist_Mem_GAPBS/Dist_Mem_GAPBS/Dist_Mem_GAPBS

module load gcc/10.1.0 openmpi/4.0.4

cd /home/zhansen/Dist_Mem_GAPBS/gapbs/src
g++ -std=c++11 -o CC cc.cc

cd /home/zhansen/Dist_Mem_GAPBS/Dist_Mem_GAPBS/Dist_Mem_GAPBS
oshc++ -std=c++11 -o CC cc.cc

#cd /home/zhansen/Dist_Mem_GAPBS/Dist_Mem_GAPBS/Dist_Mem_GAPBS
#srun -N 2 -n 4 ./CC -u 10 -n 3 -v
#mv cc_output.txt ../../cc_output.txt
#/home/zhansen/Dist_Mem_GAPBS/gapbs/src/CC -u 10 -n 3 -v
#echo "############# Trial Complete! #############\n"

#srun -N 2 -n 16 ./CC -u 15 -n 3 -v
#/home/zhansen/Dist_Mem_GAPBS/gapbs/src/CC -u 15 -n 3 -v
#echo "############# Trial Complete! #############\n"

#srun -N 2 -n 8 ./CC -g 15 -n 3 -v
#/home/zhansen/Dist_Mem_GAPBS/gapbs/src/CC -g 15 -n 3 -v
#echo "############# Trial Complete! #############\n"

#srun -N 2 -n 16 ./CC -g 15 -n 3 -v
#/home/zhansen/Dist_Mem_GAPBS/gapbs/src/CC -g 15 -n 3 -v
#echo "############# Trial Complete! #############\n"

srun -N 2 -n 4 ./CC -f /home/zhansen/Dist_Mem_GAPBS/gapbs/benchmark/graphs/raw/twitter.el -n 5 -v
mv cc_output.txt ../../cc_output.txt
/home/zhansen/Dist_Mem_GAPBS/gapbs/src/CC -f /home/zhansen/Dist_Mem_GAPBS/gapbs/benchmark/graphs/raw/twitter.el -n 5 -v
echo "############# Trial 1 Complete! #############\n"

##srun -N 1 -n 32 ./CC -f /home/zhansen/Dist_Mem_GAPBS/gapbs/benchmark/graphs/raw/twitter.el -n 1 -v
#/home/zhansen/Dist_Mem_GAPBS/gapbs/src/CC -f /home/zhansen/Dist_Mem_GAPBS/gapbs/benchmark/graphs/raw/twitter.el -n 1 -v
#echo "############# Trial 1 Complete! #############\n"
