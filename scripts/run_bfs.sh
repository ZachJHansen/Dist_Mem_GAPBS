#!/bin/bash

# run_bfs.sh
# Runs BFS kernel jobs and records results to files

# SBATCH parameters
#SBATCH --time=24:00:00
#SBATCH --nodes=16
#SBATCH --ntasks-per-node=32
#SBATCH --job-name=run_bfs
#SBATCH --output=run_bfs.%j.txt
#SBATCH --partition=standard
#SBATCH --exclusive
#SBATCH --mem=0

# Modules
module purge
module load slurm
module load gcc/10.2.0
module load PrgEnv-gnu
module unload cray-mpich
module unload cray-libsci
module load xpmem
module load udreg
module load alps
module load craype-haswell
module load htop

# Set envrionment variables for custom OMPI/UCX install
export PATH=/users/bwilliams/tt_haswell/ucx/build/bin:/users/bwilliams/tt_haswell/ompi/build/bin:$PATH
export LD_LIBRARY_PATH=/users/bwilliams/tt_haswell/ucx/build/lib:/users/bwilliams/tt_haswell/ompi/build/lib:$LD_LIBRARY_PATH
export UCX_UNIFIED_MODE=1

# Ensure results directory exists
mkdir -p $HOME/Dist_Mem_GAPBS/results/bfs

# cd to directory
cd $HOME/Dist_Mem_GAPBS

# Compile the binary
oshc++ -std=c++11 -o BFS bfs2.cc

# Run Tests
export SHMEM_SYMMETRIC_HEAP_SIZE=6656M
oshrun -display-map -np 256 -npernode 16 ./BFS -g 28 -k 64 -n 8 &> $HOME/Dist_Mem_GAPBS/results/bfs/test_1.txt
echo '===============Test 1 Complete!==============='
oshrun -display-map -np 256 -npernode 16 ./BFS -u 28 -k 64 -n 8 &> $HOME/Dist_Mem_GAPBS/results/bfs/test_2.txt
echo '===============Test 2 Complete!==============='

export SHMEM_SYMMETRIC_HEAP_SIZE=8G
oshrun -display-map -np 192 -npernode 12 ./BFS -g 28 -k 64 -n 8 &> $HOME/Dist_Mem_GAPBS/results/bfs/test_3.txt
echo '===============Test 3 Complete!==============='
oshrun -display-map -np 192 -npernode 12 ./BFS -u 28 -k 64 -n 8 &> $HOME/Dist_Mem_GAPBS/results/bfs/test_4.txt
echo '===============Test 4 Complete!==============='

export SHMEM_SYMMETRIC_HEAP_SIZE=16G
oshrun -display-map -np 64 -npernode 4 ./BFS -g 28 -k 64 -n 8 &> $HOME/Dist_Mem_GAPBS/results/bfs/test_5.txt
echo '===============Test 5 Complete!==============='
oshrun -display-map -np 64 -npernode 4 ./BFS -u 28 -k 64 -n 8 &> $HOME/Dist_Mem_GAPBS/results/bfs/test_6.txt
echo '===============Test 6 Complete!==============='

export SHMEM_SYMMETRIC_HEAP_SIZE=8G
oshrun -display-map -np 192 -npernode 12 ./BFS -g 29 -k 64 -n 8 &> $HOME/Dist_Mem_GAPBS/results/bfs/test_7.txt
echo '===============Test 7 Complete!==============='
oshrun -display-map -np 192 -npernode 12 ./BFS -u 29 -k 64 -n 8 &> $HOME/Dist_Mem_GAPBS/results/bfs/test_8.txt
echo '===============Test 8 Complete!==============='

export SHMEM_SYMMETRIC_HEAP_SIZE=18G
oshrun -display-map -np 64 -npernode 4 ./BFS -g 29 -k 64 -n 8 &> $HOME/Dist_Mem_GAPBS/results/bfs/test_9.txt
echo '===============Test 9 Complete!==============='
oshrun -display-map -np 64 -npernode 4 ./BFS -u 29 -k 64 -n 8 &> $HOME/Dist_Mem_GAPBS/results/bfs/test_10.txt
echo '===============Test 10 Complete!==============='

export SHMEM_SYMMETRIC_HEAP_SIZE=18G
oshrun -display-map -np 64 -npernode 4 ./BFS -g 30 -k 64 -n 8 &> $HOME/Dist_Mem_GAPBS/results/bfs/test_11.txt
echo '===============Test 11 Complete!==============='
oshrun -display-map -np 64 -npernode 4 ./BFS -u 30 -k 64 -n 8 &> $HOME/Dist_Mem_GAPBS/results/bfs/test_12.txt
echo '===============Test 12 Complete!==============='

export SHMEM_SYMMETRIC_HEAP_SIZE=36G
oshrun -display-map -np 32 -npernode 2 ./BFS -g 31 -k 64 -n 8 &> $HOME/Dist_Mem_GAPBS/results/bfs/test_13.txt
echo '===============Test 13 Complete!==============='
oshrun -display-map -np 32 -npernode 2 ./BFS -u 31 -k 64 -n 8 &> $HOME/Dist_Mem_GAPBS/results/bfs/test_14.txt
echo '===============Test 14 Complete!==============='
