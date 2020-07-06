// Copyright (c) 2015, The Regents of the University of California (Regents)
// See LICENSE.txt for license details

#include <algorithm>
#include <cinttypes>
#include <iostream>
#include <vector>

#include "benchmark.h"
#include "builder.h"
#include "command_line.h"
#include "graph.h"
#include "pvector.h"


/*
GAP Benchmark Suite
Kernel: Triangle Counting (TC)
Author: Scott Beamer

Will count the number of triangles (cliques of size 3)

Requires input graph:
  - to be undirected
  - no duplicate edges (or else will be counted as multiple triangles)
  - neighborhoods are sorted by vertex identifiers

Other than symmetrizing, the rest of the requirements are done by SquishCSR
during graph building.

This implementation reduces the search space by counting each triangle only
once. A naive implementation will count the same triangle six times because
each of the three vertices (u, v, w) will count it in both ways. To count
a triangle only once, this implementation only counts a triangle if u > v > w.
Once the remaining unexamined neighbors identifiers get too big, it can break
out of the loop, but this requires that the neighbors to be sorted.

Another optimization this implementation has is to relabel the vertices by
degree. This is beneficial if the average degree is high enough and if the
degree distribution is sufficiently non-uniform. To decide whether or not
to relabel the graph, we use the heuristic in WorthRelabelling.
*/


using namespace std;

size_t OrderedCount(const Graph &g, long* pSync, long* pWrk) {
  Partition<NodeID> vp(g.num_nodes());
  long* total = (long *) shmem_calloc(1, sizeof(long));
  for (NodeID u = vp.start; u < vp.end; u++) {
    for (NodeID v : g.out_neigh(u)) {
  //    printf("Pe %d | U = %d, V = %d\n", vp.pe, u, v);
      if (v > u)
        break;
      auto it = g.out_neigh(u).begin();
      for (NodeID w : g.out_neigh(v)) {
        //printf("\tV = %d, W = %d\n", v, w);
        if (w > v)
          break;
        while (*it < w) {
    //      printf("\t\tOut neigh(%d) = %d\n", u, *it);
          it++;
        }
        if (w == *it) {
      //    printf("\t\t\tIterating: u = %d, v = %d, w = %d\n", u, v, w);
          (*total)++;
        }
      }
    }
  }
//  printf("Pe %d: Total = %lu\n", vp.pe, *total);
  shmem_long_sum_to_all(total, total, 1, 0, 0, vp.npes, pWrk, pSync);                   // double should fit size_t regardless of 32 vs 64 bit system?
  return ((size_t) *total);
}

// heuristic to see if sufficently dense power-law graph                        Does this still hold for the partitioned version?
bool WorthRelabelling(const Graph &g, long *pSync, long *pWrk) {
  int64_t average_degree = g.num_edges() / g.num_nodes();
  if (average_degree < 10)
    return false;
  SourcePicker<Graph> sp(g);
  int64_t num_samples = min(int64_t(1000), g.num_nodes());
  int64_t* sample_total = (int64_t *) shmem_calloc(1, sizeof(int64_t));
  Partition<int> sample_part(num_samples);
  pvector<int64_t> samples(sample_part.max_width, true);                         // symmetric partitioned overallocated pvector
  pvector<int64_t> dest(num_samples, true);                             // symmetric overallocated pvector
  for (int64_t trial = sample_part.start; trial < sample_part.end; trial++) {
    samples[sample_part.local_pos(trial)] = g.out_degree(sp.PickNext());
    *sample_total += samples[sample_part.local_pos(trial)];
  }
  shmem_collect64(dest.begin(), samples.begin(), sample_part.end-sample_part.start, 0, 0, sample_part.npes, pSync);
  shmem_long_sum_to_all(sample_total, sample_total, 1, 0, 0, sample_part.npes, pWrk, pSync);
  sort(dest.begin(), dest.end());
  double sample_average = static_cast<double>(*sample_total) / num_samples;
  double sample_median = dest[num_samples/2];
  return sample_average / 1.3 > sample_median;
}


// uses heuristic to see if worth relabeling
size_t Hybrid(const Graph &g, long* pSync, long* pWrk) {
  if (WorthRelabelling(g, pSync, pWrk))
    return OrderedCount(Builder::RelabelByDegree(g, pSync, pWrk), pSync, pWrk);
  else
    return OrderedCount(g, pSync, pWrk);
}


void PrintTriangleStats(const Graph &g, size_t total_triangles) {
  cout << total_triangles << " triangles" << endl;
}


// Compares with simple serial implementation that uses std::set_intersection
bool TCVerifier(const Graph &g, size_t test_total) {
  if (shmem_my_pe() == 0) {
    ofstream shmem_out;
    shmem_out.open("/home/zach/projects/Dist_Mem_GAPBS/Dist_Mem_GAPBS/shmem_output.txt", ios::app);
    shmem_out << test_total << endl;
    shmem_out.close();
  }
/*
  size_t total = 0;
  vector<NodeID> intersection;
  intersection.reserve(g.num_nodes());
  for (NodeID u : g.vertices()) {
    for (NodeID v : g.out_neigh(u)) {
      auto new_end = set_intersection(g.out_neigh(u).begin(),
                                      g.out_neigh(u).end(),
                                      g.out_neigh(v).begin(),
                                      g.out_neigh(v).end(),
                                      intersection.begin());
      intersection.resize(new_end - intersection.begin());
      total += intersection.size();
    }
  }
  total = total / 6;  // each triangle was counted 6 times
  if (total != test_total)
    cout << total << " != " << test_total << endl;
  return total == test_total;*/
  return true;
}


int main(int argc, char* argv[]) {
  CLApp cli(argc, argv, "triangle count");
  if (!cli.ParseArgs())
    return -1;

  static long PRINT_LOCK = 0;

  char size_env[] = "SHMEM_SYMMETRIC_SIZE=8192M";
  putenv(size_env);

  shmem_init();

  static long pSync[SHMEM_REDUCE_SYNC_SIZE];
  static long lng_pWrk[SHMEM_REDUCE_MIN_WRKDATA_SIZE];      
  static double dbl_pWrk[SHMEM_REDUCE_MIN_WRKDATA_SIZE];                 

  for (int i = 0; i < SHMEM_REDUCE_SYNC_SIZE; i++)
    pSync[i] = SHMEM_SYNC_VALUE;
  for (int i = 0; i < SHMEM_REDUCE_MIN_WRKDATA_SIZE; i++) {
    lng_pWrk[i] = SHMEM_SYNC_VALUE;
    dbl_pWrk[i] = SHMEM_SYNC_VALUE;
  }

  {
    Builder b(cli);
//    printf("Check 1\n");
    Graph g = b.MakeGraph(pSync, lng_pWrk);
//    printf("Check 2\n");
    //g.PrintTopology();
    //Graph gee = Builder::RelabelByDegree(g, pSync, lng_pWrk);
    //gee.PrintTopology();
    if (g.directed()) {
      cout << "Input graph is directed but tc requires undirected" << endl;
      return -2;
    } 
    auto TCBound = [] (const Graph &g) { return Hybrid(g, lng_pWrk, pSync); };

    //size_t result = TCBound(gee);//Hybrid(g, pSync, lng_pWrk);
    //printf("Result: %lu\n", result);
    BenchmarkKernel(cli, g, TCBound, PrintTriangleStats, TCVerifier);
  }
  shmem_finalize();
  return 0;
}
