#!/bin/bash

# run_tc.sh
# Runs TC kernel jobs and records results to files

# SBATCH parameters
#SBATCH --time=12:00:00
#SBATCH --nodes=16
#SBATCH --ntasks-per-node=32
#SBATCH --job-name=run_tc
#SBATCH --output=run_tc.%j.txt
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
mkdir -p $HOME/Dist_Mem_GAPBS/results/tc

# cd to directory
cd $HOME/Dist_Mem_GAPBS

# Compile the binary
oshc++ -std=c++11 -o TC tc.cc

# Run Tests
export SHMEM_SYMMETRIC_HEAP_SIZE=15G
oshrun -display-map -np 256 -npernode 16 ./TC -g 28 -k 64 -n 8 &> $HOME/Dist_Mem_GAPBS/results/tc/test_1.txt
echo '===============Test 1 Complete!==============='
oshrun -display-map -np 256 -npernode 16 ./TC -u 28 -k 64 -n 8 &> $HOME/Dist_Mem_GAPBS/results/tc/test_2.txt
echo '===============Test 2 Complete!==============='

export SHMEM_SYMMETRIC_HEAP_SIZE=18G
oshrun -display-map -np 192 -npernode 12 ./TC -g 28 -k 64 -n 8 &> $HOME/Dist_Mem_GAPBS/results/tc/test_3.txt
echo '===============Test 3 Complete!==============='
oshrun -display-map -np 192 -npernode 12 ./TC -u 28 -k 64 -n 8 &> $HOME/Dist_Mem_GAPBS/results/tc/test_4.txt
echo '===============Test 4 Complete!==============='

export SHMEM_SYMMETRIC_HEAP_SIZE=28G
oshrun -display-map -np 64 -npernode 4 ./TC -g 28 -k 64 -n 8 &> $HOME/Dist_Mem_GAPBS/results/tc/test_5.txt
echo '===============Test 5 Complete!==============='
oshrun -display-map -np 64 -npernode 4 ./TC -u 28 -k 64 -n 8 &> $HOME/Dist_Mem_GAPBS/results/tc/test_6.txt
echo '===============Test 6 Complete!==============='

#export SHMEM_SYMMETRIC_HEAP_SIZE=8G
#oshrun -display-map -np 192 -npernode 12 ./TC -g 29 -k 64 -n 8 &> $HOME/Dist_Mem_GAPBS/results/tc/test_7.txt
#echo '===============Test 7 Complete!==============='
#oshrun -display-map -np 192 -npernode 12 ./TC -u 29 -k 64 -n 8 &> $HOME/Dist_Mem_GAPBS/results/tc/test_8.txt
#echo '===============Test 8 Complete!==============='

export SHMEM_SYMMETRIC_HEAP_SIZE=55G
oshrun -display-map -np 64 -npernode 4 ./TC -g 29 -k 64 -n 8 &> $HOME/Dist_Mem_GAPBS/results/tc/test_9.txt
echo '===============Test 9 Complete!==============='
oshrun -display-map -np 64 -npernode 4 ./TC -u 29 -k 64 -n 8 &> $HOME/Dist_Mem_GAPBS/results/tc/test_10.txt
echo '===============Test 10 Complete!==============='

#export SHMEM_SYMMETRIC_HEAP_SIZE=18G
#oshrun -display-map -np 64 -npernode 4 ./TC -g 30 -k 64 -n 8 &> $HOME/Dist_Mem_GAPBS/results/tc/test_11.txt
#echo '===============Test 11 Complete!==============='
#oshrun -display-map -np 64 -npernode 4 ./TC -u 30 -k 64 -n 8 &> $HOME/Dist_Mem_GAPBS/results/tc/test_12.txt
#echo '===============Test 12 Complete!==============='

#export SHMEM_SYMMETRIC_HEAP_SIZE=36G
#oshrun -display-map -np 32 -npernode 2 ./TC -g 31 -k 64 -n 8 &> $HOME/Dist_Mem_GAPBS/results/tc/test_13.txt
#echo '===============Test 13 Complete!==============='
#oshrun -display-map -np 32 -npernode 2 ./TC -u 31 -k 64 -n 8 &> $HOME/Dist_Mem_GAPBS/results/tc/test_14.txt
#echo '===============Test 14 Complete!==============='
