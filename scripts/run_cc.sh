#!/bin/bash

# run_cc.sh
# Runs CC kernel jobs and records results to files

# SBATCH parameters
#SBATCH --time=24:00:00
#SBATCH --nodes=16
#SBATCH --ntasks-per-node=32
#SBATCH --job-name=run_cc
#SBATCH --output=run_cc.%j.txt
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
mkdir -p $HOME/Dist_Mem_GAPBS/results/cc

# cd to directory
cd $HOME/Dist_Mem_GAPBS

# Compile the binary
oshc++ -std=c++11 -o CC cc.cc

# Run Tests
export SHMEM_SYMMETRIC_HEAP_SIZE=6656M # Needed 15297760892 = 14.247GB
oshrun -display-map -np 256 -npernode 16 ./CC -g 28 -k 64 -n 8 &> $HOME/Dist_Mem_GAPBS/results/cc/test_1.txt
echo '===============Test 1 Complete!==============='
oshrun -display-map -np 256 -npernode 16 ./CC -u 28 -k 64 -n 8 &> $HOME/Dist_Mem_GAPBS/results/cc/test_2.txt
echo '===============Test 2 Complete!==============='

export SHMEM_SYMMETRIC_HEAP_SIZE=8G #Needed 18714151208 = 17.428GB
oshrun -display-map -np 192 -npernode 12 ./CC -g 28 -k 64 -n 8 &> $HOME/Dist_Mem_GAPBS/results/cc/test_3.txt
echo '===============Test 3 Complete!==============='
oshrun -display-map -np 192 -npernode 12 ./CC -u 28 -k 64 -n 8 &> $HOME/Dist_Mem_GAPBS/results/cc/test_4.txt
echo '===============Test 4 Complete!==============='

export SHMEM_SYMMETRIC_HEAP_SIZE=16G #Needed 26484841272  = 24.665GB
oshrun -display-map -np 64 -npernode 4 ./CC -g 28 -k 64 -n 8 &> $HOME/Dist_Mem_GAPBS/results/cc/test_5.txt
echo '===============Test 5 Complete!==============='
oshrun -display-map -np 64 -npernode 4 ./CC -u 28 -k 64 -n 8 &> $HOME/Dist_Mem_GAPBS/results/cc/test_6.txt
echo '===============Test 6 Complete!==============='

export SHMEM_SYMMETRIC_HEAP_SIZE=8G #Needed 37431676116  = 34.860 GB
oshrun -display-map -np 192 -npernode 12 ./CC -g 29 -k 64 -n 8 &> $HOME/Dist_Mem_GAPBS/results/cc/test_7.txt
echo '===============Test 7 Complete!==============='
oshrun -display-map -np 192 -npernode 12 ./CC -u 29 -k 64 -n 8 &> $HOME/Dist_Mem_GAPBS/results/cc/test_8.txt
echo '===============Test 8 Complete!==============='

export SHMEM_SYMMETRIC_HEAP_SIZE=18G #Needed 52968840052 = 49.33 GB
oshrun -display-map -np 64 -npernode 4 ./CC -g 29 -k 64 -n 8 &> $HOME/Dist_Mem_GAPBS/results/cc/test_9.txt
echo '===============Test 9 Complete!==============='
oshrun -display-map -np 64 -npernode 4 ./CC -u 29 -k 64 -n 8 &> $HOME/Dist_Mem_GAPBS/results/cc/test_10.txt
echo '===============Test 10 Complete!==============='

export SHMEM_SYMMETRIC_HEAP_SIZE=18G #Needed 105938696540 = 98.663
oshrun -display-map -np 64 -npernode 4 ./CC -g 30 -k 64 -n 8 &> $HOME/Dist_Mem_GAPBS/results/cc/test_11.txt
echo '===============Test 11 Complete!==============='
oshrun -display-map -np 64 -npernode 4 ./CC -u 30 -k 64 -n 8 &> $HOME/Dist_Mem_GAPBS/results/cc/test_12.txt
echo '===============Test 12 Complete!==============='

export SHMEM_SYMMETRIC_HEAP_SIZE=36G #Other error
oshrun -display-map -np 32 -npernode 2 ./CC -g 31 -k 64 -n 8 &> $HOME/Dist_Mem_GAPBS/results/cc/test_13.txt
echo '===============Test 13 Complete!==============='
oshrun -display-map -np 32 -npernode 2 ./CC -u 31 -k 64 -n 8 &> $HOME/Dist_Mem_GAPBS/results/cc/test_14.txt
echo '===============Test 14 Complete!==============='
