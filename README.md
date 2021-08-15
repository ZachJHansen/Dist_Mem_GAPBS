# DMM-GAPBS: Porting the GAP Benchmark Suite to a PGAS model with OpenSHMEM and xBGAS.

This is a prototype implementation for a distributed memory model version of the GAP Benchmark Suite (https://github.com/sbeamer/gapbs). It implements a Partitioned Global Address Space (PGAS) with OpenSHMEM (http://openshmem.org/site/sites/default/site_files/OpenSHMEM-1.4.pdf) and C++11. The DMM-GAPBS is designed to process graphs that exceed the capacity of a single computer. Graphs are divided between Processing Elements (PEs) based on a naive partitioning of vertices and are built in symmetric memory. Most supporting data structures used in the kernels and in graph building are also partitioned. This reduces the space required (per PE) to store and process a graph. 

The src folder on master contains the latest stable source code. We compiled this code with GCC 10.2.0, OpenSHMEM 4.1.1 and UCX 1.10.1. Verification uses a modified version of the GAPBS (gapbs-verify). For example:

'oshc++ -std=c++11 -o BFS src/bfs.cc'           // Compiles the DMM-GAPBS implementation of Breadth-First Search

'oshrun -np 16 ./BFS -u 15 -n 5 -v'             // Executes 5 trials of BFS on a partitioned Uniform-Random graph of degree 15 with 16 PEs. Stores the result in bfs_output.txt

'g++ -std=c++11 -o BFS-V gapbs-verify/bfs.cc'   // Compiles the GAPBS BFS that verifies the result stored in bfs_output.txt

'./BFS-V -u 15 -n 5 -v'                         // Build the shared-memory graph, verify bfs_output.txt
