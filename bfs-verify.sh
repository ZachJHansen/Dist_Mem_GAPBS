#!/bin/bash
#SBATCH -J bfs-verify
#SBATCH -o %x.%j.out
#SBATCH -p nocona
#SBATCH -N 2
#SBATCH --ntasks-per-node=128
#SBATCH -t 24:00:00
#SBATCH -D /home/zhansen/Dist_Mem_GAPBS/Dist_Mem_GAPBS/Dist_Mem_GAPBS

module load gcc/10.1.0 openmpi/4.0.4

cd /home/zhansen/Dist_Mem_GAPBS/gapbs/src
g++ -std=c++11 -o BFS bfs.cc

cd /home/zhansen/Dist_Mem_GAPBS/Dist_Mem_GAPBS/Dist_Mem_GAPBS
oshc++ -std=c++11 -o BFS bfs2.cc

cd /home/zhansen/Dist_Mem_GAPBS/Dist_Mem_GAPBS/Dist_Mem_GAPBS
srun -N 2 -n 4 ./BFS -u 10 -n 3 -v
mv bfs_output.txt ../../bfs_output.txt
/home/zhansen/Dist_Mem_GAPBS/gapbs/src/BFS -u 10 -n 3 -v
echo "############# Trial Complete! #############\n"

srun -N 2 -n 16 ./BFS -u 15 -n 3 -v
mv bfs_output.txt ../../bfs_output.txt
/home/zhansen/Dist_Mem_GAPBS/gapbs/src/BFS -u 15 -n 3 -v
echo "############# Trial Complete! #############\n"

srun -N 2 -n 8 ./BFS -g 15 -n 3 -v
mv bfs_output.txt ../../bfs_output.txt
/home/zhansen/Dist_Mem_GAPBS/gapbs/src/BFS -g 15 -n 3 -v
echo "############# Trial Complete! #############\n"

srun -N 2 -n 16 ./BFS -g 15 -n 3 -v
mv bfs_output.txt ../../bfs_output.txt
/home/zhansen/Dist_Mem_GAPBS/gapbs/src/BFS -g 15 -n 3 -v
echo "############# Trial Complete! #############\n"

srun -N 2 -n 4 ./BFS -f /home/zhansen/Dist_Mem_GAPBS/gapbs/benchmark/graphs/raw/twitter.el -n 5 -v
mv bfs_output.txt ../../bfs_output.txt
/home/zhansen/Dist_Mem_GAPBS/gapbs/src/BFS -f /home/zhansen/Dist_Mem_GAPBS/gapbs/benchmark/graphs/raw/twitter.el -n 5 -v
echo "############# Trial Complete! #############\n"

srun -N 1 -n 32 ./BFS -f /home/zhansen/Dist_Mem_GAPBS/gapbs/benchmark/graphs/raw/twitter.el -n 1 -v
mv bfs_output.txt ../../bfs_output.txt
/home/zhansen/Dist_Mem_GAPBS/gapbs/src/BFS -f /home/zhansen/Dist_Mem_GAPBS/gapbs/benchmark/graphs/raw/twitter.el -n 1 -v
echo "############# Trial Complete! #############\n"
