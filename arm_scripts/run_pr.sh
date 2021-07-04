#!/bin/bash

# run_pr.sh
# Runs PR kernel jobs and records results to files

# SBATCH parameters
#SBATCH --time=12:00:00
#SBATCH --nodes=16
#SBATCH --ntasks-per-node=32
#SBATCH --job-name=run_pr
#SBATCH --output=run_pr.%j.txt
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
module load craype-arm-thunderx2
module load htop

# Set envrionment variables for custom OMPI/UCX install
export PATH=/users/bwilliams/capulin/ucx/build/bin:/users/bwilliams/capulin/ompi/build/bin:$PATH
export LD_LIBRARY_PATH=/users/bwilliams/capulin/ucx/build/lib:/users/bwilliams/capulin/ompi/build/lib:$LD_LIBRARY_PATH
export UCX_UNIFIED_MODE=1

# Ensure results directory exists
mkdir -p $HOME/Dist_Mem_GAPBS/results/pr

# cd to directory
cd $HOME/Dist_Mem_GAPBS

# Compile the binary
oshc++ -std=c++11 -o PR pr.cc

# Run Tests
export SHMEM_SYMMETRIC_HEAP_SIZE=15G
oshrun -display-map -np 256 -npernode 16 ./PR -g 28 -k 64 -n 8 &> $HOME/Dist_Mem_GAPBS/results/pr/test_1.txt
echo '===============Test 1 Complete!==============='
oshrun -display-map -np 256 -npernode 16 ./PR -u 28 -k 64 -n 8 &> $HOME/Dist_Mem_GAPBS/results/pr/test_2.txt
echo '===============Test 2 Complete!==============='

export SHMEM_SYMMETRIC_HEAP_SIZE=18G
oshrun -display-map -np 192 -npernode 12 ./PR -g 28 -k 64 -n 8 &> $HOME/Dist_Mem_GAPBS/results/pr/test_3.txt
echo '===============Test 3 Complete!==============='
oshrun -display-map -np 192 -npernode 12 ./PR -u 28 -k 64 -n 8 &> $HOME/Dist_Mem_GAPBS/results/pr/test_4.txt
echo '===============Test 4 Complete!==============='

export SHMEM_SYMMETRIC_HEAP_SIZE=28G
oshrun -display-map -np 64 -npernode 4 ./PR -g 28 -k 64 -n 8 &> $HOME/Dist_Mem_GAPBS/results/pr/test_5.txt
echo '===============Test 5 Complete!==============='
oshrun -display-map -np 64 -npernode 4 ./PR -u 28 -k 64 -n 8 &> $HOME/Dist_Mem_GAPBS/results/pr/test_6.txt
echo '===============Test 6 Complete!==============='

#export SHMEM_SYMMETRIC_HEAP_SIZE=8G
#oshrun -display-map -np 192 -npernode 12 ./PR -g 29 -k 64 -n 8 &> $HOME/Dist_Mem_GAPBS/results/pr/test_7.txt
#echo '===============Test 7 Complete!==============='
#oshrun -display-map -np 192 -npernode 12 ./PR -u 29 -k 64 -n 8 &> $HOME/Dist_Mem_GAPBS/results/pr/test_8.txt
#echo '===============Test 8 Complete!==============='

export SHMEM_SYMMETRIC_HEAP_SIZE=55G
oshrun -display-map -np 64 -npernode 4 ./PR -g 29 -k 64 -n 8 &> $HOME/Dist_Mem_GAPBS/results/pr/test_9.txt
echo '===============Test 9 Complete!==============='
oshrun -display-map -np 64 -npernode 4 ./PR -u 29 -k 64 -n 8 &> $HOME/Dist_Mem_GAPBS/results/pr/test_10.txt
echo '===============Test 10 Complete!==============='

#export SHMEM_SYMMETRIC_HEAP_SIZE=18G
#oshrun -display-map -np 64 -npernode 4 ./PR -g 30 -k 64 -n 8 &> $HOME/Dist_Mem_GAPBS/results/pr/test_11.txt
#echo '===============Test 11 Complete!==============='
#oshrun -display-map -np 64 -npernode 4 ./PR -u 30 -k 64 -n 8 &> $HOME/Dist_Mem_GAPBS/results/pr/test_12.txt
#echo '===============Test 12 Complete!==============='

#export SHMEM_SYMMETRIC_HEAP_SIZE=36G
#oshrun -display-map -np 32 -npernode 2 ./PR -g 31 -k 64 -n 8 &> $HOME/Dist_Mem_GAPBS/results/pr/test_13.txt
#echo '===============Test 13 Complete!==============='
#oshrun -display-map -np 32 -npernode 2 ./PR -u 31 -k 64 -n 8 &> $HOME/Dist_Mem_GAPBS/results/pr/test_14.txt
#echo '===============Test 14 Complete!==============='
