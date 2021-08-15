#!/bin/bash

PES="$1"
PPR="$2"
TOTAL_PES="$3"
DEGREE="$4"
TEST="BC"
HOSTFILE="output-$TEST-$PES-$PPR.txt"

cd /home/zhansen/Dist_Mem_GAPBS
source ./set_env.sh
oshc++ -std=c++11 -o BC bc.cc 
oshrun -display-map -np $TOTAL_PES --npernode $PPR ./BC -g $DEGREE -n 3 -v
../gapbs/src/BC -g $DEGREE -n 3 -v
