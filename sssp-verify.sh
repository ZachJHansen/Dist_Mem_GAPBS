#!/bin/bash
#SBATCH -J sssp-verify
#SBATCH -o %x.%j.out
#SBATCH -p nocona
#SBATCH -N 8
#SBATCH --ntasks-per-node=128
#SBATCH -t 120:00:00
#SBATCH -D /home/zhansen/Dist_Mem_GAPBS/Dist_Mem_GAPBS/Dist_Mem_GAPBS

module load gcc/10.1.0 openmpi/4.0.4

cd /home/zhansen/Dist_Mem_GAPBS/gapbs/src
g++ -std=c++11 -o SSSP sssp.cc

cd /home/zhansen/Dist_Mem_GAPBS/Dist_Mem_GAPBS/Dist_Mem_GAPBS
oshc++ -std=c++11 -o SSSP new_sssp.cc

srun -N 2 -n 4 ./SSSP -u 20 -d 10 -n 16 -v
mv sssp_output.txt ../../sssp_output.txt
/home/zhansen/Dist_Mem_GAPBS/gapbs/src/SSSP -u 20 -d 10 -n 16 -v
echo '############# Trial 1 Complete! #############'

srun -N 2 -n 4 ./SSSP -g 20 -d 10 -n 16 -v
mv sssp_output.txt ../../sssp_output.txt
/home/zhansen/Dist_Mem_GAPBS/gapbs/src/SSSP -g 20 -d 10 -n 16 -v
echo '############# Trial 1 Complete! #############'

srun -N 2 -n 4 ./SSSP -u 20 -d 1000 -n 16 -v
mv sssp_output.txt ../../sssp_output.txt
/home/zhansen/Dist_Mem_GAPBS/gapbs/src/SSSP -u 20 -d 1000 -n 16 -v
echo '############# Trial 2 Complete! #############'

srun -N 2 -n 4 ./SSSP -g 20 -d 1000 -n 16 -v
mv sssp_output.txt ../../sssp_output.txt
/home/zhansen/Dist_Mem_GAPBS/gapbs/src/SSSP -g 20 -d 1000 -n 16 -v
echo '############# Trial 2 Complete! #############'

srun -N 2 -n 4 ./SSSP -u 20 -d 10000 -n 16 -v
mv sssp_output.txt ../../sssp_output.txt
/home/zhansen/Dist_Mem_GAPBS/gapbs/src/SSSP -u 20 -d 10000 -n 16 -v
echo '############# Trial 3 Complete! #############'

srun -N 2 -n 4 ./SSSP -g 20 -d 10000 -n 16 -v
mv sssp_output.txt ../../sssp_output.txt
/home/zhansen/Dist_Mem_GAPBS/gapbs/src/SSSP -g 20 -d 10000 -n 16 -v
echo '############# Trial 3 Complete! #############'

srun -N 2 -n 4 ./SSSP -u 20 -d 1000000 -n 16 -v
mv sssp_output.txt ../../sssp_output.txt
/home/zhansen/Dist_Mem_GAPBS/gapbs/src/SSSP -u 20 -d 1000000 -n 16 -v
echo '############# Trial 4 Complete! #############'

srun -N 2 -n 4 ./SSSP -g 20 -d 1000000 -n 16 -v
mv sssp_output.txt ../../sssp_output.txt
/home/zhansen/Dist_Mem_GAPBS/gapbs/src/SSSP -g 20 -d 1000000 -n 16 -v
echo '############# Trial 4 Complete! #############'

srun -N 2 -n 4 ./SSSP -u 30 -d 10 -n 16 -v
mv sssp_output.txt ../../sssp_output.txt
/home/zhansen/Dist_Mem_GAPBS/gapbs/src/SSSP -u 30 -d 10 -n 16 -v
echo '############# Trial 5 Complete! #############'

srun -N 2 -n 4 ./SSSP -g 30 -d 10 -n 16 -v
mv sssp_output.txt ../../sssp_output.txt
/home/zhansen/Dist_Mem_GAPBS/gapbs/src/SSSP -g 30 -d 10 -n 16 -v
echo '############# Trial 5 Complete! #############'

srun -N 2 -n 4 ./SSSP -u 30 -d 1000 -n 16 -v
mv sssp_output.txt ../../sssp_output.txt
/home/zhansen/Dist_Mem_GAPBS/gapbs/src/SSSP -u 30 -d 1000 -n 16 -v
echo '############# Trial 6 Complete! #############'

srun -N 2 -n 4 ./SSSP -g 30 -d 1000 -n 16 -v
mv sssp_output.txt ../../sssp_output.txt
/home/zhansen/Dist_Mem_GAPBS/gapbs/src/SSSP -g 30 -d 1000 -n 16 -v
echo '############# Trial 6 Complete! #############'

srun -N 2 -n 4 ./SSSP -u 30 -d 10000 -n 16 -v
mv sssp_output.txt ../../sssp_output.txt
/home/zhansen/Dist_Mem_GAPBS/gapbs/src/SSSP -u 30 -d 10000 -n 16 -v
echo '############# Trial 7 Complete! #############'

srun -N 2 -n 4 ./SSSP -g 30 -d 10000 -n 16 -v
mv sssp_output.txt ../../sssp_output.txt
/home/zhansen/Dist_Mem_GAPBS/gapbs/src/SSSP -g 30 -d 10000 -n 16 -v
echo '############# Trial 7 Complete! #############'

srun -N 2 -n 4 ./SSSP -u 30 -d 1000000 -n 16 -v
mv sssp_output.txt ../../sssp_output.txt
/home/zhansen/Dist_Mem_GAPBS/gapbs/src/SSSP -u 30 -d 1000000 -n 16 -v
echo '############# Trial 8 Complete! #############'

srun -N 2 -n 4 ./SSSP -g 30 -d 1000000 -n 16 -v
mv sssp_output.txt ../../sssp_output.txt
/home/zhansen/Dist_Mem_GAPBS/gapbs/src/SSSP -g 30 -d 1000000 -n 16 -v
echo '############# Trial 8 Complete! #############'

srun -N 2 -n 4 ./SSSP -f /home/zhansen/Dist_Mem_GAPBS/gapbs/benchmark/graphs/raw/twitter.el -d 10 -n 16 -v
mv sssp_output.txt ../../sssp_output.txt
/home/zhansen/Dist_Mem_GAPBS/gapbs/src/SSSP -f /home/zhansen/Dist_Mem_GAPBS/gapbs/benchmark/graphs/raw/twitter.el -d 10 -n 16 -v
echo '############# Trial 9 Complete! #############'

srun -N 2 -n 4 ./SSSP -f /home/zhansen/Dist_Mem_GAPBS/gapbs/benchmark/graphs/raw/twitter.el -d 1000 -n 16 -v
mv sssp_output.txt ../../sssp_output.txt
/home/zhansen/Dist_Mem_GAPBS/gapbs/src/SSSP -f /home/zhansen/Dist_Mem_GAPBS/gapbs/benchmark/graphs/raw/twitter.el -d 1000 -n 16 -v
echo '############# Trial 10 Complete! #############'

srun -N 2 -n 4 ./SSSP -f /home/zhansen/Dist_Mem_GAPBS/gapbs/benchmark/graphs/raw/twitter.el -d 10000 -n 16 -v
mv sssp_output.txt ../../sssp_output.txt
/home/zhansen/Dist_Mem_GAPBS/gapbs/src/SSSP -f /home/zhansen/Dist_Mem_GAPBS/gapbs/benchmark/graphs/raw/twitter.el -d 10000 -n 16 -v
echo '############# Trial 11 Complete! #############'

srun -N 2 -n 4 ./SSSP -f /home/zhansen/Dist_Mem_GAPBS/gapbs/benchmark/graphs/raw/twitter.el -d 1000000 -n 16 -v
mv sssp_output.txt ../../sssp_output.txt
/home/zhansen/Dist_Mem_GAPBS/gapbs/src/SSSP -f /home/zhansen/Dist_Mem_GAPBS/gapbs/benchmark/graphs/raw/twitter.el -d 1000000 -n 16 -v
echo '############# Trial 12 Complete! #############'

srun -N 2 -n 4 ./SSSP -f /home/zhansen/Dist_Mem_GAPBS/gapbs/benchmark/graphs/raw/USA-road-d.USA.gr -d 10 -n 16 -v
mv sssp_output.txt ../../sssp_output.txt
/home/zhansen/Dist_Mem_GAPBS/gapbs/src/SSSP -f /home/zhansen/Dist_Mem_GAPBS/gapbs/benchmark/graphs/raw/USA-road-d.USA.gr -d 10 -n 16 -v
echo '############# Trial 13 Complete! #############'

srun -N 2 -n 4 ./SSSP -f /home/zhansen/Dist_Mem_GAPBS/gapbs/benchmark/graphs/raw/USA-road-d.USA.gr -d 1000 -n 16 -v
mv sssp_output.txt ../../sssp_output.txt
/home/zhansen/Dist_Mem_GAPBS/gapbs/src/SSSP -f /home/zhansen/Dist_Mem_GAPBS/gapbs/benchmark/graphs/raw/USA-road-d.USA.gr -d 1000 -n 16 -v
echo '############# Trial 14 Complete! #############'

srun -N 2 -n 4 ./SSSP -f /home/zhansen/Dist_Mem_GAPBS/gapbs/benchmark/graphs/raw/USA-road-d.USA.gr -d 10000 -n 16 -v
mv sssp_output.txt ../../sssp_output.txt
/home/zhansen/Dist_Mem_GAPBS/gapbs/src/SSSP -f /home/zhansen/Dist_Mem_GAPBS/gapbs/benchmark/graphs/raw/USA-road-d.USA.gr -d 10000 -n 16 -v
echo '############# Trial 15 Complete! #############'

srun -N 2 -n 4 ./SSSP -f /home/zhansen/Dist_Mem_GAPBS/gapbs/benchmark/graphs/raw/USA-road-d.USA.gr -d 1000000 -n 16 -v
mv sssp_output.txt ../../sssp_output.txt
/home/zhansen/Dist_Mem_GAPBS/gapbs/src/SSSP -f /home/zhansen/Dist_Mem_GAPBS/gapbs/benchmark/graphs/raw/USA-road-d.USA.gr -d 1000000 -n 16 -v
echo '############# Trial 16 Complete! #############'

srun -N 2 -n 16 ./SSSP -u 20 -d 10 -n 16 -v
mv sssp_output.txt ../../sssp_output.txt
/home/zhansen/Dist_Mem_GAPBS/gapbs/src/SSSP -u 20 -d 10 -n 16 -v
echo '############# Trial 17 Complete! #############'

srun -N 2 -n 16 ./SSSP -g 20 -d 10 -n 16 -v
mv sssp_output.txt ../../sssp_output.txt
/home/zhansen/Dist_Mem_GAPBS/gapbs/src/SSSP -g 20 -d 10 -n 16 -v
echo '############# Trial 17 Complete! #############'

srun -N 2 -n 16 ./SSSP -u 20 -d 1000 -n 16 -v
mv sssp_output.txt ../../sssp_output.txt
/home/zhansen/Dist_Mem_GAPBS/gapbs/src/SSSP -u 20 -d 1000 -n 16 -v
echo '############# Trial 18 Complete! #############'

srun -N 2 -n 16 ./SSSP -g 20 -d 1000 -n 16 -v
mv sssp_output.txt ../../sssp_output.txt
/home/zhansen/Dist_Mem_GAPBS/gapbs/src/SSSP -g 20 -d 1000 -n 16 -v
echo '############# Trial 18 Complete! #############'

srun -N 2 -n 16 ./SSSP -u 20 -d 10000 -n 16 -v
mv sssp_output.txt ../../sssp_output.txt
/home/zhansen/Dist_Mem_GAPBS/gapbs/src/SSSP -u 20 -d 10000 -n 16 -v
echo '############# Trial 19 Complete! #############'

srun -N 2 -n 16 ./SSSP -g 20 -d 10000 -n 16 -v
mv sssp_output.txt ../../sssp_output.txt
/home/zhansen/Dist_Mem_GAPBS/gapbs/src/SSSP -g 20 -d 10000 -n 16 -v
echo '############# Trial 19 Complete! #############'

srun -N 2 -n 16 ./SSSP -u 20 -d 1000000 -n 16 -v
mv sssp_output.txt ../../sssp_output.txt
/home/zhansen/Dist_Mem_GAPBS/gapbs/src/SSSP -u 20 -d 1000000 -n 16 -v
echo '############# Trial 20 Complete! #############'

srun -N 2 -n 16 ./SSSP -g 20 -d 1000000 -n 16 -v
mv sssp_output.txt ../../sssp_output.txt
/home/zhansen/Dist_Mem_GAPBS/gapbs/src/SSSP -g 20 -d 1000000 -n 16 -v
echo '############# Trial 20 Complete! #############'

srun -N 2 -n 16 ./SSSP -u 30 -d 10 -n 16 -v
mv sssp_output.txt ../../sssp_output.txt
/home/zhansen/Dist_Mem_GAPBS/gapbs/src/SSSP -u 30 -d 10 -n 16 -v
echo '############# Trial 21 Complete! #############'

srun -N 2 -n 16 ./SSSP -g 30 -d 10 -n 16 -v
mv sssp_output.txt ../../sssp_output.txt
/home/zhansen/Dist_Mem_GAPBS/gapbs/src/SSSP -g 30 -d 10 -n 16 -v
echo '############# Trial 21 Complete! #############'

srun -N 2 -n 16 ./SSSP -u 30 -d 1000 -n 16 -v
mv sssp_output.txt ../../sssp_output.txt
/home/zhansen/Dist_Mem_GAPBS/gapbs/src/SSSP -u 30 -d 1000 -n 16 -v
echo '############# Trial 22 Complete! #############'

srun -N 2 -n 16 ./SSSP -g 30 -d 1000 -n 16 -v
mv sssp_output.txt ../../sssp_output.txt
/home/zhansen/Dist_Mem_GAPBS/gapbs/src/SSSP -g 30 -d 1000 -n 16 -v
echo '############# Trial 22 Complete! #############'

srun -N 2 -n 16 ./SSSP -u 30 -d 10000 -n 16 -v
mv sssp_output.txt ../../sssp_output.txt
/home/zhansen/Dist_Mem_GAPBS/gapbs/src/SSSP -u 30 -d 10000 -n 16 -v
echo '############# Trial 23 Complete! #############'

srun -N 2 -n 16 ./SSSP -g 30 -d 10000 -n 16 -v
mv sssp_output.txt ../../sssp_output.txt
/home/zhansen/Dist_Mem_GAPBS/gapbs/src/SSSP -g 30 -d 10000 -n 16 -v
echo '############# Trial 23 Complete! #############'

srun -N 2 -n 16 ./SSSP -u 30 -d 1000000 -n 16 -v
mv sssp_output.txt ../../sssp_output.txt
/home/zhansen/Dist_Mem_GAPBS/gapbs/src/SSSP -u 30 -d 1000000 -n 16 -v
echo '############# Trial 24 Complete! #############'

srun -N 2 -n 16 ./SSSP -g 30 -d 1000000 -n 16 -v
mv sssp_output.txt ../../sssp_output.txt
/home/zhansen/Dist_Mem_GAPBS/gapbs/src/SSSP -g 30 -d 1000000 -n 16 -v
echo '############# Trial 24 Complete! #############'

srun -N 2 -n 16 ./SSSP -f /home/zhansen/Dist_Mem_GAPBS/gapbs/benchmark/graphs/raw/twitter.el -d 10 -n 16 -v
mv sssp_output.txt ../../sssp_output.txt
/home/zhansen/Dist_Mem_GAPBS/gapbs/src/SSSP -f /home/zhansen/Dist_Mem_GAPBS/gapbs/benchmark/graphs/raw/twitter.el -d 10 -n 16 -v
echo '############# Trial 25 Complete! #############'

srun -N 2 -n 16 ./SSSP -f /home/zhansen/Dist_Mem_GAPBS/gapbs/benchmark/graphs/raw/twitter.el -d 1000 -n 16 -v
mv sssp_output.txt ../../sssp_output.txt
/home/zhansen/Dist_Mem_GAPBS/gapbs/src/SSSP -f /home/zhansen/Dist_Mem_GAPBS/gapbs/benchmark/graphs/raw/twitter.el -d 1000 -n 16 -v
echo '############# Trial 26 Complete! #############'

srun -N 2 -n 16 ./SSSP -f /home/zhansen/Dist_Mem_GAPBS/gapbs/benchmark/graphs/raw/twitter.el -d 10000 -n 16 -v
mv sssp_output.txt ../../sssp_output.txt
/home/zhansen/Dist_Mem_GAPBS/gapbs/src/SSSP -f /home/zhansen/Dist_Mem_GAPBS/gapbs/benchmark/graphs/raw/twitter.el -d 10000 -n 16 -v
echo '############# Trial 27 Complete! #############'

srun -N 2 -n 16 ./SSSP -f /home/zhansen/Dist_Mem_GAPBS/gapbs/benchmark/graphs/raw/twitter.el -d 1000000 -n 16 -v
mv sssp_output.txt ../../sssp_output.txt
/home/zhansen/Dist_Mem_GAPBS/gapbs/src/SSSP -f /home/zhansen/Dist_Mem_GAPBS/gapbs/benchmark/graphs/raw/twitter.el -d 1000000 -n 16 -v
echo '############# Trial 28 Complete! #############'

srun -N 2 -n 16 ./SSSP -f /home/zhansen/Dist_Mem_GAPBS/gapbs/benchmark/graphs/raw/USA-road-d.USA.gr -d 10 -n 16 -v
mv sssp_output.txt ../../sssp_output.txt
/home/zhansen/Dist_Mem_GAPBS/gapbs/src/SSSP -f /home/zhansen/Dist_Mem_GAPBS/gapbs/benchmark/graphs/raw/USA-road-d.USA.gr -d 10 -n 16 -v
echo '############# Trial 29 Complete! #############'

srun -N 2 -n 16 ./SSSP -f /home/zhansen/Dist_Mem_GAPBS/gapbs/benchmark/graphs/raw/USA-road-d.USA.gr -d 1000 -n 16 -v
mv sssp_output.txt ../../sssp_output.txt
/home/zhansen/Dist_Mem_GAPBS/gapbs/src/SSSP -f /home/zhansen/Dist_Mem_GAPBS/gapbs/benchmark/graphs/raw/USA-road-d.USA.gr -d 1000 -n 16 -v
echo '############# Trial 30 Complete! #############'

srun -N 2 -n 16 ./SSSP -f /home/zhansen/Dist_Mem_GAPBS/gapbs/benchmark/graphs/raw/USA-road-d.USA.gr -d 10000 -n 16 -v
mv sssp_output.txt ../../sssp_output.txt
/home/zhansen/Dist_Mem_GAPBS/gapbs/src/SSSP -f /home/zhansen/Dist_Mem_GAPBS/gapbs/benchmark/graphs/raw/USA-road-d.USA.gr -d 10000 -n 16 -v
echo '############# Trial 31 Complete! #############'

srun -N 2 -n 16 ./SSSP -f /home/zhansen/Dist_Mem_GAPBS/gapbs/benchmark/graphs/raw/USA-road-d.USA.gr -d 1000000 -n 16 -v
mv sssp_output.txt ../../sssp_output.txt
/home/zhansen/Dist_Mem_GAPBS/gapbs/src/SSSP -f /home/zhansen/Dist_Mem_GAPBS/gapbs/benchmark/graphs/raw/USA-road-d.USA.gr -d 1000000 -n 16 -v
echo '############# Trial 32 Complete! #############'

srun -N 4 -n 4 ./SSSP -u 20 -d 10 -n 16 -v
mv sssp_output.txt ../../sssp_output.txt
/home/zhansen/Dist_Mem_GAPBS/gapbs/src/SSSP -u 20 -d 10 -n 16 -v
echo '############# Trial 33 Complete! #############'

srun -N 4 -n 4 ./SSSP -g 20 -d 10 -n 16 -v
mv sssp_output.txt ../../sssp_output.txt
/home/zhansen/Dist_Mem_GAPBS/gapbs/src/SSSP -g 20 -d 10 -n 16 -v
echo '############# Trial 33 Complete! #############'

srun -N 4 -n 4 ./SSSP -u 20 -d 1000 -n 16 -v
mv sssp_output.txt ../../sssp_output.txt
/home/zhansen/Dist_Mem_GAPBS/gapbs/src/SSSP -u 20 -d 1000 -n 16 -v
echo '############# Trial 34 Complete! #############'

srun -N 4 -n 4 ./SSSP -g 20 -d 1000 -n 16 -v
mv sssp_output.txt ../../sssp_output.txt
/home/zhansen/Dist_Mem_GAPBS/gapbs/src/SSSP -g 20 -d 1000 -n 16 -v
echo '############# Trial 34 Complete! #############'

srun -N 4 -n 4 ./SSSP -u 20 -d 10000 -n 16 -v
mv sssp_output.txt ../../sssp_output.txt
/home/zhansen/Dist_Mem_GAPBS/gapbs/src/SSSP -u 20 -d 10000 -n 16 -v
echo '############# Trial 35 Complete! #############'

srun -N 4 -n 4 ./SSSP -g 20 -d 10000 -n 16 -v
mv sssp_output.txt ../../sssp_output.txt
/home/zhansen/Dist_Mem_GAPBS/gapbs/src/SSSP -g 20 -d 10000 -n 16 -v
echo '############# Trial 35 Complete! #############'

srun -N 4 -n 4 ./SSSP -u 20 -d 1000000 -n 16 -v
mv sssp_output.txt ../../sssp_output.txt
/home/zhansen/Dist_Mem_GAPBS/gapbs/src/SSSP -u 20 -d 1000000 -n 16 -v
echo '############# Trial 36 Complete! #############'

srun -N 4 -n 4 ./SSSP -g 20 -d 1000000 -n 16 -v
mv sssp_output.txt ../../sssp_output.txt
/home/zhansen/Dist_Mem_GAPBS/gapbs/src/SSSP -g 20 -d 1000000 -n 16 -v
echo '############# Trial 36 Complete! #############'

srun -N 4 -n 4 ./SSSP -u 30 -d 10 -n 16 -v
mv sssp_output.txt ../../sssp_output.txt
/home/zhansen/Dist_Mem_GAPBS/gapbs/src/SSSP -u 30 -d 10 -n 16 -v
echo '############# Trial 37 Complete! #############'

srun -N 4 -n 4 ./SSSP -g 30 -d 10 -n 16 -v
mv sssp_output.txt ../../sssp_output.txt
/home/zhansen/Dist_Mem_GAPBS/gapbs/src/SSSP -g 30 -d 10 -n 16 -v
echo '############# Trial 37 Complete! #############'

srun -N 4 -n 4 ./SSSP -u 30 -d 1000 -n 16 -v
mv sssp_output.txt ../../sssp_output.txt
/home/zhansen/Dist_Mem_GAPBS/gapbs/src/SSSP -u 30 -d 1000 -n 16 -v
echo '############# Trial 38 Complete! #############'

srun -N 4 -n 4 ./SSSP -g 30 -d 1000 -n 16 -v
mv sssp_output.txt ../../sssp_output.txt
/home/zhansen/Dist_Mem_GAPBS/gapbs/src/SSSP -g 30 -d 1000 -n 16 -v
echo '############# Trial 38 Complete! #############'

srun -N 4 -n 4 ./SSSP -u 30 -d 10000 -n 16 -v
mv sssp_output.txt ../../sssp_output.txt
/home/zhansen/Dist_Mem_GAPBS/gapbs/src/SSSP -u 30 -d 10000 -n 16 -v
echo '############# Trial 39 Complete! #############'

srun -N 4 -n 4 ./SSSP -g 30 -d 10000 -n 16 -v
mv sssp_output.txt ../../sssp_output.txt
/home/zhansen/Dist_Mem_GAPBS/gapbs/src/SSSP -g 30 -d 10000 -n 16 -v
echo '############# Trial 39 Complete! #############'

srun -N 4 -n 4 ./SSSP -u 30 -d 1000000 -n 16 -v
mv sssp_output.txt ../../sssp_output.txt
/home/zhansen/Dist_Mem_GAPBS/gapbs/src/SSSP -u 30 -d 1000000 -n 16 -v
echo '############# Trial 40 Complete! #############'

srun -N 4 -n 4 ./SSSP -g 30 -d 1000000 -n 16 -v
mv sssp_output.txt ../../sssp_output.txt
/home/zhansen/Dist_Mem_GAPBS/gapbs/src/SSSP -g 30 -d 1000000 -n 16 -v
echo '############# Trial 40 Complete! #############'

srun -N 4 -n 4 ./SSSP -f /home/zhansen/Dist_Mem_GAPBS/gapbs/benchmark/graphs/raw/twitter.el -d 10 -n 16 -v
mv sssp_output.txt ../../sssp_output.txt
/home/zhansen/Dist_Mem_GAPBS/gapbs/src/SSSP -f /home/zhansen/Dist_Mem_GAPBS/gapbs/benchmark/graphs/raw/twitter.el -d 10 -n 16 -v
echo '############# Trial 41 Complete! #############'

srun -N 4 -n 4 ./SSSP -f /home/zhansen/Dist_Mem_GAPBS/gapbs/benchmark/graphs/raw/twitter.el -d 1000 -n 16 -v
mv sssp_output.txt ../../sssp_output.txt
/home/zhansen/Dist_Mem_GAPBS/gapbs/src/SSSP -f /home/zhansen/Dist_Mem_GAPBS/gapbs/benchmark/graphs/raw/twitter.el -d 1000 -n 16 -v
echo '############# Trial 42 Complete! #############'

srun -N 4 -n 4 ./SSSP -f /home/zhansen/Dist_Mem_GAPBS/gapbs/benchmark/graphs/raw/twitter.el -d 10000 -n 16 -v
mv sssp_output.txt ../../sssp_output.txt
/home/zhansen/Dist_Mem_GAPBS/gapbs/src/SSSP -f /home/zhansen/Dist_Mem_GAPBS/gapbs/benchmark/graphs/raw/twitter.el -d 10000 -n 16 -v
echo '############# Trial 43 Complete! #############'

srun -N 4 -n 4 ./SSSP -f /home/zhansen/Dist_Mem_GAPBS/gapbs/benchmark/graphs/raw/twitter.el -d 1000000 -n 16 -v
mv sssp_output.txt ../../sssp_output.txt
/home/zhansen/Dist_Mem_GAPBS/gapbs/src/SSSP -f /home/zhansen/Dist_Mem_GAPBS/gapbs/benchmark/graphs/raw/twitter.el -d 1000000 -n 16 -v
echo '############# Trial 44 Complete! #############'

srun -N 4 -n 4 ./SSSP -f /home/zhansen/Dist_Mem_GAPBS/gapbs/benchmark/graphs/raw/USA-road-d.USA.gr -d 10 -n 16 -v
mv sssp_output.txt ../../sssp_output.txt
/home/zhansen/Dist_Mem_GAPBS/gapbs/src/SSSP -f /home/zhansen/Dist_Mem_GAPBS/gapbs/benchmark/graphs/raw/USA-road-d.USA.gr -d 10 -n 16 -v
echo '############# Trial 45 Complete! #############'

srun -N 4 -n 4 ./SSSP -f /home/zhansen/Dist_Mem_GAPBS/gapbs/benchmark/graphs/raw/USA-road-d.USA.gr -d 1000 -n 16 -v
mv sssp_output.txt ../../sssp_output.txt
/home/zhansen/Dist_Mem_GAPBS/gapbs/src/SSSP -f /home/zhansen/Dist_Mem_GAPBS/gapbs/benchmark/graphs/raw/USA-road-d.USA.gr -d 1000 -n 16 -v
echo '############# Trial 46 Complete! #############'

srun -N 4 -n 4 ./SSSP -f /home/zhansen/Dist_Mem_GAPBS/gapbs/benchmark/graphs/raw/USA-road-d.USA.gr -d 10000 -n 16 -v
mv sssp_output.txt ../../sssp_output.txt
/home/zhansen/Dist_Mem_GAPBS/gapbs/src/SSSP -f /home/zhansen/Dist_Mem_GAPBS/gapbs/benchmark/graphs/raw/USA-road-d.USA.gr -d 10000 -n 16 -v
echo '############# Trial 47 Complete! #############'

srun -N 4 -n 4 ./SSSP -f /home/zhansen/Dist_Mem_GAPBS/gapbs/benchmark/graphs/raw/USA-road-d.USA.gr -d 1000000 -n 16 -v
mv sssp_output.txt ../../sssp_output.txt
/home/zhansen/Dist_Mem_GAPBS/gapbs/src/SSSP -f /home/zhansen/Dist_Mem_GAPBS/gapbs/benchmark/graphs/raw/USA-road-d.USA.gr -d 1000000 -n 16 -v
echo '############# Trial 48 Complete! #############'

srun -N 4 -n 16 ./SSSP -u 20 -d 10 -n 16 -v
mv sssp_output.txt ../../sssp_output.txt
/home/zhansen/Dist_Mem_GAPBS/gapbs/src/SSSP -u 20 -d 10 -n 16 -v
echo '############# Trial 49 Complete! #############'

srun -N 4 -n 16 ./SSSP -g 20 -d 10 -n 16 -v
mv sssp_output.txt ../../sssp_output.txt
/home/zhansen/Dist_Mem_GAPBS/gapbs/src/SSSP -g 20 -d 10 -n 16 -v
echo '############# Trial 49 Complete! #############'

srun -N 4 -n 16 ./SSSP -u 20 -d 1000 -n 16 -v
mv sssp_output.txt ../../sssp_output.txt
/home/zhansen/Dist_Mem_GAPBS/gapbs/src/SSSP -u 20 -d 1000 -n 16 -v
echo '############# Trial 50 Complete! #############'

srun -N 4 -n 16 ./SSSP -g 20 -d 1000 -n 16 -v
mv sssp_output.txt ../../sssp_output.txt
/home/zhansen/Dist_Mem_GAPBS/gapbs/src/SSSP -g 20 -d 1000 -n 16 -v
echo '############# Trial 50 Complete! #############'

srun -N 4 -n 16 ./SSSP -u 20 -d 10000 -n 16 -v
mv sssp_output.txt ../../sssp_output.txt
/home/zhansen/Dist_Mem_GAPBS/gapbs/src/SSSP -u 20 -d 10000 -n 16 -v
echo '############# Trial 51 Complete! #############'

srun -N 4 -n 16 ./SSSP -g 20 -d 10000 -n 16 -v
mv sssp_output.txt ../../sssp_output.txt
/home/zhansen/Dist_Mem_GAPBS/gapbs/src/SSSP -g 20 -d 10000 -n 16 -v
echo '############# Trial 51 Complete! #############'

srun -N 4 -n 16 ./SSSP -u 20 -d 1000000 -n 16 -v
mv sssp_output.txt ../../sssp_output.txt
/home/zhansen/Dist_Mem_GAPBS/gapbs/src/SSSP -u 20 -d 1000000 -n 16 -v
echo '############# Trial 52 Complete! #############'

srun -N 4 -n 16 ./SSSP -g 20 -d 1000000 -n 16 -v
mv sssp_output.txt ../../sssp_output.txt
/home/zhansen/Dist_Mem_GAPBS/gapbs/src/SSSP -g 20 -d 1000000 -n 16 -v
echo '############# Trial 52 Complete! #############'

srun -N 4 -n 16 ./SSSP -u 30 -d 10 -n 16 -v
mv sssp_output.txt ../../sssp_output.txt
/home/zhansen/Dist_Mem_GAPBS/gapbs/src/SSSP -u 30 -d 10 -n 16 -v
echo '############# Trial 53 Complete! #############'

srun -N 4 -n 16 ./SSSP -g 30 -d 10 -n 16 -v
mv sssp_output.txt ../../sssp_output.txt
/home/zhansen/Dist_Mem_GAPBS/gapbs/src/SSSP -g 30 -d 10 -n 16 -v
echo '############# Trial 53 Complete! #############'

srun -N 4 -n 16 ./SSSP -u 30 -d 1000 -n 16 -v
mv sssp_output.txt ../../sssp_output.txt
/home/zhansen/Dist_Mem_GAPBS/gapbs/src/SSSP -u 30 -d 1000 -n 16 -v
echo '############# Trial 54 Complete! #############'

srun -N 4 -n 16 ./SSSP -g 30 -d 1000 -n 16 -v
mv sssp_output.txt ../../sssp_output.txt
/home/zhansen/Dist_Mem_GAPBS/gapbs/src/SSSP -g 30 -d 1000 -n 16 -v
echo '############# Trial 54 Complete! #############'

srun -N 4 -n 16 ./SSSP -u 30 -d 10000 -n 16 -v
mv sssp_output.txt ../../sssp_output.txt
/home/zhansen/Dist_Mem_GAPBS/gapbs/src/SSSP -u 30 -d 10000 -n 16 -v
echo '############# Trial 55 Complete! #############'

srun -N 4 -n 16 ./SSSP -g 30 -d 10000 -n 16 -v
mv sssp_output.txt ../../sssp_output.txt
/home/zhansen/Dist_Mem_GAPBS/gapbs/src/SSSP -g 30 -d 10000 -n 16 -v
echo '############# Trial 55 Complete! #############'

srun -N 4 -n 16 ./SSSP -u 30 -d 1000000 -n 16 -v
mv sssp_output.txt ../../sssp_output.txt
/home/zhansen/Dist_Mem_GAPBS/gapbs/src/SSSP -u 30 -d 1000000 -n 16 -v
echo '############# Trial 56 Complete! #############'

srun -N 4 -n 16 ./SSSP -g 30 -d 1000000 -n 16 -v
mv sssp_output.txt ../../sssp_output.txt
/home/zhansen/Dist_Mem_GAPBS/gapbs/src/SSSP -g 30 -d 1000000 -n 16 -v
echo '############# Trial 56 Complete! #############'

srun -N 4 -n 16 ./SSSP -f /home/zhansen/Dist_Mem_GAPBS/gapbs/benchmark/graphs/raw/twitter.el -d 10 -n 16 -v
mv sssp_output.txt ../../sssp_output.txt
/home/zhansen/Dist_Mem_GAPBS/gapbs/src/SSSP -f /home/zhansen/Dist_Mem_GAPBS/gapbs/benchmark/graphs/raw/twitter.el -d 10 -n 16 -v
echo '############# Trial 57 Complete! #############'

srun -N 4 -n 16 ./SSSP -f /home/zhansen/Dist_Mem_GAPBS/gapbs/benchmark/graphs/raw/twitter.el -d 1000 -n 16 -v
mv sssp_output.txt ../../sssp_output.txt
/home/zhansen/Dist_Mem_GAPBS/gapbs/src/SSSP -f /home/zhansen/Dist_Mem_GAPBS/gapbs/benchmark/graphs/raw/twitter.el -d 1000 -n 16 -v
echo '############# Trial 58 Complete! #############'

srun -N 4 -n 16 ./SSSP -f /home/zhansen/Dist_Mem_GAPBS/gapbs/benchmark/graphs/raw/twitter.el -d 10000 -n 16 -v
mv sssp_output.txt ../../sssp_output.txt
/home/zhansen/Dist_Mem_GAPBS/gapbs/src/SSSP -f /home/zhansen/Dist_Mem_GAPBS/gapbs/benchmark/graphs/raw/twitter.el -d 10000 -n 16 -v
echo '############# Trial 59 Complete! #############'

srun -N 4 -n 16 ./SSSP -f /home/zhansen/Dist_Mem_GAPBS/gapbs/benchmark/graphs/raw/twitter.el -d 1000000 -n 16 -v
mv sssp_output.txt ../../sssp_output.txt
/home/zhansen/Dist_Mem_GAPBS/gapbs/src/SSSP -f /home/zhansen/Dist_Mem_GAPBS/gapbs/benchmark/graphs/raw/twitter.el -d 1000000 -n 16 -v
echo '############# Trial 60 Complete! #############'

srun -N 4 -n 16 ./SSSP -f /home/zhansen/Dist_Mem_GAPBS/gapbs/benchmark/graphs/raw/USA-road-d.USA.gr -d 10 -n 16 -v
mv sssp_output.txt ../../sssp_output.txt
/home/zhansen/Dist_Mem_GAPBS/gapbs/src/SSSP -f /home/zhansen/Dist_Mem_GAPBS/gapbs/benchmark/graphs/raw/USA-road-d.USA.gr -d 10 -n 16 -v
echo '############# Trial 61 Complete! #############'

srun -N 4 -n 16 ./SSSP -f /home/zhansen/Dist_Mem_GAPBS/gapbs/benchmark/graphs/raw/USA-road-d.USA.gr -d 1000 -n 16 -v
mv sssp_output.txt ../../sssp_output.txt
/home/zhansen/Dist_Mem_GAPBS/gapbs/src/SSSP -f /home/zhansen/Dist_Mem_GAPBS/gapbs/benchmark/graphs/raw/USA-road-d.USA.gr -d 1000 -n 16 -v
echo '############# Trial 62 Complete! #############'

srun -N 4 -n 16 ./SSSP -f /home/zhansen/Dist_Mem_GAPBS/gapbs/benchmark/graphs/raw/USA-road-d.USA.gr -d 10000 -n 16 -v
mv sssp_output.txt ../../sssp_output.txt
/home/zhansen/Dist_Mem_GAPBS/gapbs/src/SSSP -f /home/zhansen/Dist_Mem_GAPBS/gapbs/benchmark/graphs/raw/USA-road-d.USA.gr -d 10000 -n 16 -v
echo '############# Trial 63 Complete! #############'

srun -N 4 -n 16 ./SSSP -f /home/zhansen/Dist_Mem_GAPBS/gapbs/benchmark/graphs/raw/USA-road-d.USA.gr -d 1000000 -n 16 -v
mv sssp_output.txt ../../sssp_output.txt
/home/zhansen/Dist_Mem_GAPBS/gapbs/src/SSSP -f /home/zhansen/Dist_Mem_GAPBS/gapbs/benchmark/graphs/raw/USA-road-d.USA.gr -d 1000000 -n 16 -v
echo '############# Trial 64 Complete! #############'

srun -N 8 -n 4 ./SSSP -u 20 -d 10 -n 16 -v
mv sssp_output.txt ../../sssp_output.txt
/home/zhansen/Dist_Mem_GAPBS/gapbs/src/SSSP -u 20 -d 10 -n 16 -v
echo '############# Trial 65 Complete! #############'

srun -N 8 -n 4 ./SSSP -g 20 -d 10 -n 16 -v
mv sssp_output.txt ../../sssp_output.txt
/home/zhansen/Dist_Mem_GAPBS/gapbs/src/SSSP -g 20 -d 10 -n 16 -v
echo '############# Trial 65 Complete! #############'

srun -N 8 -n 4 ./SSSP -u 20 -d 1000 -n 16 -v
mv sssp_output.txt ../../sssp_output.txt
/home/zhansen/Dist_Mem_GAPBS/gapbs/src/SSSP -u 20 -d 1000 -n 16 -v
echo '############# Trial 66 Complete! #############'

srun -N 8 -n 4 ./SSSP -g 20 -d 1000 -n 16 -v
mv sssp_output.txt ../../sssp_output.txt
/home/zhansen/Dist_Mem_GAPBS/gapbs/src/SSSP -g 20 -d 1000 -n 16 -v
echo '############# Trial 66 Complete! #############'

srun -N 8 -n 4 ./SSSP -u 20 -d 10000 -n 16 -v
mv sssp_output.txt ../../sssp_output.txt
/home/zhansen/Dist_Mem_GAPBS/gapbs/src/SSSP -u 20 -d 10000 -n 16 -v
echo '############# Trial 67 Complete! #############'

srun -N 8 -n 4 ./SSSP -g 20 -d 10000 -n 16 -v
mv sssp_output.txt ../../sssp_output.txt
/home/zhansen/Dist_Mem_GAPBS/gapbs/src/SSSP -g 20 -d 10000 -n 16 -v
echo '############# Trial 67 Complete! #############'

srun -N 8 -n 4 ./SSSP -u 20 -d 1000000 -n 16 -v
mv sssp_output.txt ../../sssp_output.txt
/home/zhansen/Dist_Mem_GAPBS/gapbs/src/SSSP -u 20 -d 1000000 -n 16 -v
echo '############# Trial 68 Complete! #############'

srun -N 8 -n 4 ./SSSP -g 20 -d 1000000 -n 16 -v
mv sssp_output.txt ../../sssp_output.txt
/home/zhansen/Dist_Mem_GAPBS/gapbs/src/SSSP -g 20 -d 1000000 -n 16 -v
echo '############# Trial 68 Complete! #############'

srun -N 8 -n 4 ./SSSP -u 30 -d 10 -n 16 -v
mv sssp_output.txt ../../sssp_output.txt
/home/zhansen/Dist_Mem_GAPBS/gapbs/src/SSSP -u 30 -d 10 -n 16 -v
echo '############# Trial 69 Complete! #############'

srun -N 8 -n 4 ./SSSP -g 30 -d 10 -n 16 -v
mv sssp_output.txt ../../sssp_output.txt
/home/zhansen/Dist_Mem_GAPBS/gapbs/src/SSSP -g 30 -d 10 -n 16 -v
echo '############# Trial 69 Complete! #############'

srun -N 8 -n 4 ./SSSP -u 30 -d 1000 -n 16 -v
mv sssp_output.txt ../../sssp_output.txt
/home/zhansen/Dist_Mem_GAPBS/gapbs/src/SSSP -u 30 -d 1000 -n 16 -v
echo '############# Trial 70 Complete! #############'

srun -N 8 -n 4 ./SSSP -g 30 -d 1000 -n 16 -v
mv sssp_output.txt ../../sssp_output.txt
/home/zhansen/Dist_Mem_GAPBS/gapbs/src/SSSP -g 30 -d 1000 -n 16 -v
echo '############# Trial 70 Complete! #############'

srun -N 8 -n 4 ./SSSP -u 30 -d 10000 -n 16 -v
mv sssp_output.txt ../../sssp_output.txt
/home/zhansen/Dist_Mem_GAPBS/gapbs/src/SSSP -u 30 -d 10000 -n 16 -v
echo '############# Trial 71 Complete! #############'

srun -N 8 -n 4 ./SSSP -g 30 -d 10000 -n 16 -v
mv sssp_output.txt ../../sssp_output.txt
/home/zhansen/Dist_Mem_GAPBS/gapbs/src/SSSP -g 30 -d 10000 -n 16 -v
echo '############# Trial 71 Complete! #############'

srun -N 8 -n 4 ./SSSP -u 30 -d 1000000 -n 16 -v
mv sssp_output.txt ../../sssp_output.txt
/home/zhansen/Dist_Mem_GAPBS/gapbs/src/SSSP -u 30 -d 1000000 -n 16 -v
echo '############# Trial 72 Complete! #############'

srun -N 8 -n 4 ./SSSP -g 30 -d 1000000 -n 16 -v
mv sssp_output.txt ../../sssp_output.txt
/home/zhansen/Dist_Mem_GAPBS/gapbs/src/SSSP -g 30 -d 1000000 -n 16 -v
echo '############# Trial 72 Complete! #############'

srun -N 8 -n 4 ./SSSP -f /home/zhansen/Dist_Mem_GAPBS/gapbs/benchmark/graphs/raw/twitter.el -d 10 -n 16 -v
mv sssp_output.txt ../../sssp_output.txt
/home/zhansen/Dist_Mem_GAPBS/gapbs/src/SSSP -f /home/zhansen/Dist_Mem_GAPBS/gapbs/benchmark/graphs/raw/twitter.el -d 10 -n 16 -v
echo '############# Trial 73 Complete! #############'

srun -N 8 -n 4 ./SSSP -f /home/zhansen/Dist_Mem_GAPBS/gapbs/benchmark/graphs/raw/twitter.el -d 1000 -n 16 -v
mv sssp_output.txt ../../sssp_output.txt
/home/zhansen/Dist_Mem_GAPBS/gapbs/src/SSSP -f /home/zhansen/Dist_Mem_GAPBS/gapbs/benchmark/graphs/raw/twitter.el -d 1000 -n 16 -v
echo '############# Trial 74 Complete! #############'

srun -N 8 -n 4 ./SSSP -f /home/zhansen/Dist_Mem_GAPBS/gapbs/benchmark/graphs/raw/twitter.el -d 10000 -n 16 -v
mv sssp_output.txt ../../sssp_output.txt
/home/zhansen/Dist_Mem_GAPBS/gapbs/src/SSSP -f /home/zhansen/Dist_Mem_GAPBS/gapbs/benchmark/graphs/raw/twitter.el -d 10000 -n 16 -v
echo '############# Trial 75 Complete! #############'

srun -N 8 -n 4 ./SSSP -f /home/zhansen/Dist_Mem_GAPBS/gapbs/benchmark/graphs/raw/twitter.el -d 1000000 -n 16 -v
mv sssp_output.txt ../../sssp_output.txt
/home/zhansen/Dist_Mem_GAPBS/gapbs/src/SSSP -f /home/zhansen/Dist_Mem_GAPBS/gapbs/benchmark/graphs/raw/twitter.el -d 1000000 -n 16 -v
echo '############# Trial 76 Complete! #############'

srun -N 8 -n 4 ./SSSP -f /home/zhansen/Dist_Mem_GAPBS/gapbs/benchmark/graphs/raw/USA-road-d.USA.gr -d 10 -n 16 -v
mv sssp_output.txt ../../sssp_output.txt
/home/zhansen/Dist_Mem_GAPBS/gapbs/src/SSSP -f /home/zhansen/Dist_Mem_GAPBS/gapbs/benchmark/graphs/raw/USA-road-d.USA.gr -d 10 -n 16 -v
echo '############# Trial 77 Complete! #############'

srun -N 8 -n 4 ./SSSP -f /home/zhansen/Dist_Mem_GAPBS/gapbs/benchmark/graphs/raw/USA-road-d.USA.gr -d 1000 -n 16 -v
mv sssp_output.txt ../../sssp_output.txt
/home/zhansen/Dist_Mem_GAPBS/gapbs/src/SSSP -f /home/zhansen/Dist_Mem_GAPBS/gapbs/benchmark/graphs/raw/USA-road-d.USA.gr -d 1000 -n 16 -v
echo '############# Trial 78 Complete! #############'

srun -N 8 -n 4 ./SSSP -f /home/zhansen/Dist_Mem_GAPBS/gapbs/benchmark/graphs/raw/USA-road-d.USA.gr -d 10000 -n 16 -v
mv sssp_output.txt ../../sssp_output.txt
/home/zhansen/Dist_Mem_GAPBS/gapbs/src/SSSP -f /home/zhansen/Dist_Mem_GAPBS/gapbs/benchmark/graphs/raw/USA-road-d.USA.gr -d 10000 -n 16 -v
echo '############# Trial 79 Complete! #############'

srun -N 8 -n 4 ./SSSP -f /home/zhansen/Dist_Mem_GAPBS/gapbs/benchmark/graphs/raw/USA-road-d.USA.gr -d 1000000 -n 16 -v
mv sssp_output.txt ../../sssp_output.txt
/home/zhansen/Dist_Mem_GAPBS/gapbs/src/SSSP -f /home/zhansen/Dist_Mem_GAPBS/gapbs/benchmark/graphs/raw/USA-road-d.USA.gr -d 1000000 -n 16 -v
echo '############# Trial 80 Complete! #############'

srun -N 8 -n 16 ./SSSP -u 20 -d 10 -n 16 -v
mv sssp_output.txt ../../sssp_output.txt
/home/zhansen/Dist_Mem_GAPBS/gapbs/src/SSSP -u 20 -d 10 -n 16 -v
echo '############# Trial 81 Complete! #############'

srun -N 8 -n 16 ./SSSP -g 20 -d 10 -n 16 -v
mv sssp_output.txt ../../sssp_output.txt
/home/zhansen/Dist_Mem_GAPBS/gapbs/src/SSSP -g 20 -d 10 -n 16 -v
echo '############# Trial 81 Complete! #############'

srun -N 8 -n 16 ./SSSP -u 20 -d 1000 -n 16 -v
mv sssp_output.txt ../../sssp_output.txt
/home/zhansen/Dist_Mem_GAPBS/gapbs/src/SSSP -u 20 -d 1000 -n 16 -v
echo '############# Trial 82 Complete! #############'

srun -N 8 -n 16 ./SSSP -g 20 -d 1000 -n 16 -v
mv sssp_output.txt ../../sssp_output.txt
/home/zhansen/Dist_Mem_GAPBS/gapbs/src/SSSP -g 20 -d 1000 -n 16 -v
echo '############# Trial 82 Complete! #############'

srun -N 8 -n 16 ./SSSP -u 20 -d 10000 -n 16 -v
mv sssp_output.txt ../../sssp_output.txt
/home/zhansen/Dist_Mem_GAPBS/gapbs/src/SSSP -u 20 -d 10000 -n 16 -v
echo '############# Trial 83 Complete! #############'

srun -N 8 -n 16 ./SSSP -g 20 -d 10000 -n 16 -v
mv sssp_output.txt ../../sssp_output.txt
/home/zhansen/Dist_Mem_GAPBS/gapbs/src/SSSP -g 20 -d 10000 -n 16 -v
echo '############# Trial 83 Complete! #############'

srun -N 8 -n 16 ./SSSP -u 20 -d 1000000 -n 16 -v
mv sssp_output.txt ../../sssp_output.txt
/home/zhansen/Dist_Mem_GAPBS/gapbs/src/SSSP -u 20 -d 1000000 -n 16 -v
echo '############# Trial 84 Complete! #############'

srun -N 8 -n 16 ./SSSP -g 20 -d 1000000 -n 16 -v
mv sssp_output.txt ../../sssp_output.txt
/home/zhansen/Dist_Mem_GAPBS/gapbs/src/SSSP -g 20 -d 1000000 -n 16 -v
echo '############# Trial 84 Complete! #############'

srun -N 8 -n 16 ./SSSP -u 30 -d 10 -n 16 -v
mv sssp_output.txt ../../sssp_output.txt
/home/zhansen/Dist_Mem_GAPBS/gapbs/src/SSSP -u 30 -d 10 -n 16 -v
echo '############# Trial 85 Complete! #############'

srun -N 8 -n 16 ./SSSP -g 30 -d 10 -n 16 -v
mv sssp_output.txt ../../sssp_output.txt
/home/zhansen/Dist_Mem_GAPBS/gapbs/src/SSSP -g 30 -d 10 -n 16 -v
echo '############# Trial 85 Complete! #############'

srun -N 8 -n 16 ./SSSP -u 30 -d 1000 -n 16 -v
mv sssp_output.txt ../../sssp_output.txt
/home/zhansen/Dist_Mem_GAPBS/gapbs/src/SSSP -u 30 -d 1000 -n 16 -v
echo '############# Trial 86 Complete! #############'

srun -N 8 -n 16 ./SSSP -g 30 -d 1000 -n 16 -v
mv sssp_output.txt ../../sssp_output.txt
/home/zhansen/Dist_Mem_GAPBS/gapbs/src/SSSP -g 30 -d 1000 -n 16 -v
echo '############# Trial 86 Complete! #############'

srun -N 8 -n 16 ./SSSP -u 30 -d 10000 -n 16 -v
mv sssp_output.txt ../../sssp_output.txt
/home/zhansen/Dist_Mem_GAPBS/gapbs/src/SSSP -u 30 -d 10000 -n 16 -v
echo '############# Trial 87 Complete! #############'

srun -N 8 -n 16 ./SSSP -g 30 -d 10000 -n 16 -v
mv sssp_output.txt ../../sssp_output.txt
/home/zhansen/Dist_Mem_GAPBS/gapbs/src/SSSP -g 30 -d 10000 -n 16 -v
echo '############# Trial 87 Complete! #############'

srun -N 8 -n 16 ./SSSP -u 30 -d 1000000 -n 16 -v
mv sssp_output.txt ../../sssp_output.txt
/home/zhansen/Dist_Mem_GAPBS/gapbs/src/SSSP -u 30 -d 1000000 -n 16 -v
echo '############# Trial 88 Complete! #############'

srun -N 8 -n 16 ./SSSP -g 30 -d 1000000 -n 16 -v
mv sssp_output.txt ../../sssp_output.txt
/home/zhansen/Dist_Mem_GAPBS/gapbs/src/SSSP -g 30 -d 1000000 -n 16 -v
echo '############# Trial 88 Complete! #############'

srun -N 8 -n 16 ./SSSP -f /home/zhansen/Dist_Mem_GAPBS/gapbs/benchmark/graphs/raw/twitter.el -d 10 -n 16 -v
mv sssp_output.txt ../../sssp_output.txt
/home/zhansen/Dist_Mem_GAPBS/gapbs/src/SSSP -f /home/zhansen/Dist_Mem_GAPBS/gapbs/benchmark/graphs/raw/twitter.el -d 10 -n 16 -v
echo '############# Trial 89 Complete! #############'

srun -N 8 -n 16 ./SSSP -f /home/zhansen/Dist_Mem_GAPBS/gapbs/benchmark/graphs/raw/twitter.el -d 1000 -n 16 -v
mv sssp_output.txt ../../sssp_output.txt
/home/zhansen/Dist_Mem_GAPBS/gapbs/src/SSSP -f /home/zhansen/Dist_Mem_GAPBS/gapbs/benchmark/graphs/raw/twitter.el -d 1000 -n 16 -v
echo '############# Trial 90 Complete! #############'

srun -N 8 -n 16 ./SSSP -f /home/zhansen/Dist_Mem_GAPBS/gapbs/benchmark/graphs/raw/twitter.el -d 10000 -n 16 -v
mv sssp_output.txt ../../sssp_output.txt
/home/zhansen/Dist_Mem_GAPBS/gapbs/src/SSSP -f /home/zhansen/Dist_Mem_GAPBS/gapbs/benchmark/graphs/raw/twitter.el -d 10000 -n 16 -v
echo '############# Trial 91 Complete! #############'

srun -N 8 -n 16 ./SSSP -f /home/zhansen/Dist_Mem_GAPBS/gapbs/benchmark/graphs/raw/twitter.el -d 1000000 -n 16 -v
mv sssp_output.txt ../../sssp_output.txt
/home/zhansen/Dist_Mem_GAPBS/gapbs/src/SSSP -f /home/zhansen/Dist_Mem_GAPBS/gapbs/benchmark/graphs/raw/twitter.el -d 1000000 -n 16 -v
echo '############# Trial 92 Complete! #############'

srun -N 8 -n 16 ./SSSP -f /home/zhansen/Dist_Mem_GAPBS/gapbs/benchmark/graphs/raw/USA-road-d.USA.gr -d 10 -n 16 -v
mv sssp_output.txt ../../sssp_output.txt
/home/zhansen/Dist_Mem_GAPBS/gapbs/src/SSSP -f /home/zhansen/Dist_Mem_GAPBS/gapbs/benchmark/graphs/raw/USA-road-d.USA.gr -d 10 -n 16 -v
echo '############# Trial 93 Complete! #############'

srun -N 8 -n 16 ./SSSP -f /home/zhansen/Dist_Mem_GAPBS/gapbs/benchmark/graphs/raw/USA-road-d.USA.gr -d 1000 -n 16 -v
mv sssp_output.txt ../../sssp_output.txt
/home/zhansen/Dist_Mem_GAPBS/gapbs/src/SSSP -f /home/zhansen/Dist_Mem_GAPBS/gapbs/benchmark/graphs/raw/USA-road-d.USA.gr -d 1000 -n 16 -v
echo '############# Trial 94 Complete! #############'

srun -N 8 -n 16 ./SSSP -f /home/zhansen/Dist_Mem_GAPBS/gapbs/benchmark/graphs/raw/USA-road-d.USA.gr -d 10000 -n 16 -v
mv sssp_output.txt ../../sssp_output.txt
/home/zhansen/Dist_Mem_GAPBS/gapbs/src/SSSP -f /home/zhansen/Dist_Mem_GAPBS/gapbs/benchmark/graphs/raw/USA-road-d.USA.gr -d 10000 -n 16 -v
echo '############# Trial 95 Complete! #############'

srun -N 8 -n 16 ./SSSP -f /home/zhansen/Dist_Mem_GAPBS/gapbs/benchmark/graphs/raw/USA-road-d.USA.gr -d 1000000 -n 16 -v
mv sssp_output.txt ../../sssp_output.txt
/home/zhansen/Dist_Mem_GAPBS/gapbs/src/SSSP -f /home/zhansen/Dist_Mem_GAPBS/gapbs/benchmark/graphs/raw/USA-road-d.USA.gr -d 1000000 -n 16 -v
echo '############# Trial 96 Complete! #############'
