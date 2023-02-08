#!/bin/sh
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=4
#SBATCH --mem=96gb
#SBATCH --time=16:00:00
#SBATCH --job-name=testing
#SBATCH --error=testing.%J.err
#SBATCH --output=testing.%J.out

module load gcc/10.1 openmpi/4.0
oshc++ -std=c++11 -o BC src/bc.cc
oshc++ -std=c++11 -o CC src/cc.cc
oshc++ -std=c++11 -o BFS src/bfs.cc
oshc++ -std=c++11 -o TC src/tc.cc
oshc++ -std=c++11 -o PR src/pr.cc
oshc++ -std=c++11 -o SSSP src/sssp.cc

python test_harness.py BC all

python test_harness.py CC all

python test_harness.py BFS all

python test_harness.py TC all

python test_harness.py PR all

python test_harness.py SSSP all
