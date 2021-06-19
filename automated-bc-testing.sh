#!/bin/bash

TIME="4:00:00"
PART="defq"
PERNODE="3"
NODES="2"
DEGREE="10"
cd /home/zhansen/Dist_Mem_GAPBS
for PPN in $PERNODE;do
    for NODE in $NODES;do
        for DEG in $DEGREE;do
            TASKS=$( bc -l <<<"$PPN*$NODE" )
            echo "BC Test: SUBMITTING TOTAL_TASKS=$TASKS; NODES=$NODE; TASKS PER NODE=$PPN; GRAPH DEGREE=$DEG"
            sbatch -o BC-SYNTHETIC-$NODE-$PPN.out --partition=$PART --time=$TIME --ntasks=$TASKS --nodes=$NODE --ntasks-per-node=$PPN --exclusive BC-TEST.sh $NODE $PPN $TASKS $DEG
        done;
    done;
done;
