#!/bin/bash

PES="$1"
PPR="$2"
TOTAL_PES="$3"
DEGREE="$4"
TEST="PR"
HOSTFILE="output-$TEST-$PES-$PPR.txt"

cd /home/zhansen/gapbs/src
g++ -std=c++11 -o PR pr.cc

cd /home/zhansen/Dist_Mem_GAPBS
source ./set_env.sh
oshc++ -std=c++11 -o PR pr.cc 
oshrun -display-map -np $TOTAL_PES --npernode $PPR ./PR -g $DEGREE -n 3 -v
../gapbs/src/PR -g $DEGREE -n 3 -v
