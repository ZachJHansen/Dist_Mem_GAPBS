#!/bin/sh

oshc++ -std=c++11 -o BC bc.cc
python test_harness.py BC all

oshc++ -std=c++11 -o CC cc.cc
python test_harness.py CC all

oshc++ -std=c++11 -o BFS bfs.cc
python test_harness.py BFS all

oshc++ -std=c++11 -o TC tc.cc
python test_harness.py TC all

oshc++ -std=c++11 -o PR pr.cc
python test_harness.py PR all

oshc++ -std=c++11 -o SSSP original_sssp.cc
python test_harness.py SSSP all
