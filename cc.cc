// Copyright (c) 2018, The Hebrew University of Jerusalem (HUJI, A. Barak)
// See LICENSE.txt for license details

#include <algorithm>
#include <cinttypes>
#include <iostream>
#include <unordered_map>
#include <vector>

#include "benchmark.h"
#include "bitmap.h"
#include "builder.h"
#include "command_line.h"
#include "graph.h"
#include "pvector.h"
#include "timer.h"


/*
GAP Benchmark Suite
Kernel: Connected Components (CC)
Authors: Michael Sutton, Scott Beamer

Will return comp array labelling each vertex with a connected component ID

This CC implementation makes use of the Afforest subgraph sampling algorithm [1],
which restructures and extends the Shiloach-Vishkin algorithm [2].

[1] Michael Sutton, Tal Ben-Nun, and Amnon Barak. "Optimizing Parallel 
    Graph Connectivity Computation via Subgraph Sampling" Symposium on 
    Parallel and Distributed Processing, IPDPS 2018.

[2] Yossi Shiloach and Uzi Vishkin. "An o(logn) parallel connectivity algorithm"
    Journal of Algorithms, 3(1):57–67, 1982.
*/




// so [1] is saying that min-label propogation is a good algorithm for distributed memory envs
using namespace std;


// Place nodes u and v in same component of lower component ID  
void Link(NodeID u, NodeID v, pvector<int>& comp, Partition<NodeID> vp) {   
  int temp, comparator;
  int p1 = shmem_int_g(comp.begin()+vp.local_pos(u), vp.recv(u));
  int p2 = shmem_int_g(comp.begin()+vp.local_pos(v), vp.recv(v));
  while (p1 != p2) {
    int high = p1 > p2 ? p1 : p2;
    int low = p1 + (p2 - high);
    int p_high = shmem_int_g(comp.begin()+vp.local_pos(high), vp.recv(high));
    // Was already 'low' or succeeded in writing 'low'
    //if ((p_high == low) ||
        //(p_high == high && shmem_int_atomic_compare_swap(comp.begin()+vp.local_pos(high), high, low, vp.recv(high))))
    if ((p_high == low) || ((p_high == high) && (shmem_int_atomic_compare_swap(comp.begin()+vp.local_pos(high), high, low, vp.recv(high)) == high)))
      break;
    /*if (p_high == low) {
      break;
    } else {
      int prev = shmem_int_atomic_compare_swap(comp.begin()+vp.local_pos(high), high, low, vp.recv(high));
      if ((p_high == high) && (prev == high))
        break;
    }*/
    temp = shmem_int_g(comp.begin()+vp.local_pos(high), vp.recv(high));
    p1 = shmem_int_g(comp.begin()+vp.local_pos(temp), vp.recv(temp));
    p2 = shmem_int_g(comp.begin()+vp.local_pos(low), vp.recv(low));
  }
}


// Reduce depth of tree for each component to 1 by crawling up parents
void Compress(const Graph &g, pvector<int>& comp, Partition<NodeID> cp) {
  /*#pragma omp parallel for schedule(dynamic, 16384)
  for (NodeID n = 0; n < g.num_nodes(); n++) {
    while (comp[n] != comp[comp[n]]) {
      comp[n] = comp[comp[n]];
    }
  }*/
  //int pop_of_n, own_of_p, l_of_p, loc_n, p_of_n;
  for (NodeID n = cp.start; n < cp.end; n++) {
    while(comp[cp.local_pos(n)] != shmem_int_g(comp.begin() + cp.local_pos(comp[cp.local_pos(n)]), cp.recv(comp[cp.local_pos(n)]))) {
      comp[cp.local_pos(n)] = shmem_int_g(comp.begin() + cp.local_pos(comp[cp.local_pos(n)]), cp.recv(comp[cp.local_pos(n)]));
    }

/*
    loc_n = cp.local_pos(n);                    // local position of n
    while (true) {
      p_of_n = comp[loc_n];                     // parent of n (guaranteed to be local since n is within bounds of cp)
      l_of_p = cp.local_pos(p_of_n);            // local position of parent
      own_of_p = cp.recv(p_of_n);               // PE to whom parent is assigned
      pop_of_n = shmem_int_g(comp.begin()+l_of_p, own_of_p);    // parent of parent of n
      comp[loc_n] = pop_of_n;
      //printf("PE %d | N = %d, p = %d, own(p) = %d, l(p) = %d, pop(n) = %d\n", cp.pe, n, p_of_n, own_of_p, l_of_p, pop_of_n);
      if (comp[loc_n] == pop_of_n)
        break;
    }*/

  }
}

// so if i partition & parrallelize this so each node contributes 1024/npes samples, we lose some randomness
// also what if npes > 1024?
// maybe its best for the first pe that reaches this function to handle it alone, since work is small
void SampleFrequentElement(const pvector<int>& comp, int* most_freq, Partition<NodeID> vp,
                             int64_t num_samples = 1024) {
  std::unordered_map<NodeID, int> sample_counts(32);                    // 32 pairs, key = NodeID, val = long
  using kvp_type = std::unordered_map<NodeID, int>::value_type;         // kvp_type is a pair (2 tuple) of type <NodeID, int>
  // Sample elements from 'comp'
  std::mt19937 gen;
  std::uniform_int_distribution<NodeID> distribution(0, vp.N - 1);       // ints from 0 to num_nodes-1 with equal prob
  for (NodeID i = 0; i < num_samples; i++) {
    NodeID n = distribution(gen);
    sample_counts[shmem_int_g(comp.begin()+vp.local_pos(n), vp.recv(n))]++;
    //sample_counts[comp[n]]++;           // increment value of key val pair with key comp[n]
  }
  // Find most frequent element in samples (estimate of most frequent overall)
  auto most_frequent = std::max_element(
    sample_counts.begin(), sample_counts.end(),
    [](const kvp_type& a, const kvp_type& b) { return a.second < b.second; });
  float frac_of_graph = static_cast<float>(most_frequent->second) / num_samples;
  std::cout
    << "Skipping largest intermediate component (ID: " << most_frequent->first
    << ", approx. " << static_cast<int>(frac_of_graph * 100)
    << "% of the graph)" << std::endl;
  for (int i = 0; i < vp.npes; i++)
    if (i != vp.pe)
      shmem_int_p(most_freq, most_frequent->first, i);
  //return most_frequent->first;
}


pvector<int> Afforest(const Graph &g, long* pSync, int32_t neighbor_rounds = 2) {
  Partition<NodeID> vp(g.num_nodes());
  pvector<int> comp(vp.max_width, true);                       // symmetric partitioned (unsynched) array

  // Initialize each node to a single-node self-pointing tree
  for (NodeID n = vp.start; n < vp.end; n++)
    comp[vp.local_pos(n)] = n;                                  // label should be global "name" (index) of node, even if position is local
  // Process a sparse sampled subgraph first for approximating components.
  // Sample by processing a fixed number of neighbors for each node (see paper)
  shmem_barrier_all();
  for (int r = 0; r < neighbor_rounds; ++r) {
    for (NodeID u = vp.start; u < vp.end; u++) {
      for (NodeID v : g.out_neigh(u, r)) {
        // Link at most one time if neighbor available at offset r
        Link(u, v, comp, vp);
        break;
      }
    }
    shmem_barrier_all();                // necessary?
    Compress(g, comp, vp);
    shmem_barrier_all();
  }
  // Sample 'comp' to find the most frequent element -- due to prior
  // compression, this value represents the largest intermediate component
  int* c = (int*) shmem_malloc(sizeof(int));
  if (vp.pe == 0)
    SampleFrequentElement(comp, c, vp);
  shmem_barrier_all();

  // Final 'link' phase over remaining edges (excluding largest component)
  if (!g.directed()) {
    //#pragma omp parallel for schedule(dynamic, 16384)
    for (NodeID u = vp.start; u < vp.end; u++) {
      // Skip processing nodes in the largest component
      if (shmem_int_g(comp.begin()+vp.local_pos(u), vp.recv(u)) == *c)
        continue;
      // Skip over part of neighborhood (determined by neighbor_rounds)
      for (NodeID v : g.out_neigh(u, neighbor_rounds)) {
        Link(u, v, comp, vp);
      }
    }
  } else {
    //#pragma omp parallel for schedule(dynamic, 16384)
    for (NodeID u = vp.start; u < vp.end; u++) {
      if (shmem_int_g(comp.begin()+vp.local_pos(u), vp.recv(u)) == *c)
        continue;
      for (NodeID v : g.out_neigh(u, neighbor_rounds)) {
        Link(u, v, comp, vp);
      }
      // To support directed graphs, process reverse graph completely
      for (NodeID v : g.in_neigh(u)) {
        Link(u, v, comp, vp);
      }
    }
  }
  // Finally, 'compress' for final convergence
  shmem_barrier_all();                          // necessary?
  Compress(g, comp, vp);
  /*shmem_barrier_all();
  pvector<int> labels = comp.combine(g.num_nodes(), pSync);
  shmem_barrier_all();
  return labels;*/
  return comp;
}


void PrintCompStats(const Graph &g, const pvector<int> &comp) {
  cout << endl;
  unordered_map<int, int> count;
  for (int comp_i : comp)
    count[comp_i] += 1;
  int k = 5;
  vector<pair<int, int>> count_vector;
  count_vector.reserve(count.size());
  for (auto kvp : count)
    count_vector.push_back(kvp);
  vector<pair<int, int>> top_k = TopK(count_vector, k);
  k = min(k, static_cast<int>(top_k.size()));
  cout << k << " biggest clusters" << endl;
  for (auto kvp : top_k)
    cout << kvp.second << ":" << kvp.first << endl;
  cout << "There are " << count.size() << " components" << endl;
}


// Verifies CC result by performing a BFS from a vertex in each component
// - Asserts search does not reach a vertex with a different component label
// - If the graph is directed, it performs the search as if it was undirected
// - Asserts every vertex is visited (degree-0 vertex should have own label)
bool CCVerifier(const Graph &g, const pvector<int> &comp) {
  Partition<NodeID> vp(g.num_nodes());
  int* PRINTER = (int *) shmem_malloc(sizeof(int));
  *PRINTER = 0;
  shmem_barrier_all();
  shmem_int_wait_until(PRINTER, SHMEM_CMP_EQ, vp.pe);           // wait until previous PE puts your pe # in PRINTER
  /*printf("comp begin: %p\n", (void*) comp.begin());
  if (vp.pe == 0) {
    ofstream shmem_out;
    shmem_out.open("/home/zach/projects/Dist_Mem_GAPBS/Dist_Mem_GAPBS/cc_output.txt", ios::app);
    for (long lable : comp)
      shmem_out << lable << endl;
    */
  ofstream shmem_out;
  shmem_out.open("/home/zhansen/Dist_Mem_GAPBS/Dist_Mem_GAPBS/Dist_Mem_GAPBS/cc_output.txt", ios::app);
  for (NodeID n = vp.start; n < vp.end; n++) 
    shmem_out << comp[vp.local_pos(n)] << endl;
  shmem_out.close();
  //}
  if (!(vp.pe == vp.npes-1))
    shmem_int_p(PRINTER, vp.pe+1, vp.pe+1);             // who's next?
  return true;
}


int main(int argc, char* argv[]) {
  CLApp cli(argc, argv, "connected-components-afforest");
  if (!cli.ParseArgs())
    return -1;   

  char size_env[] = "SMA_SYMMETRIC_SIZE=16G";
  putenv(size_env);

  shmem_init();

  static long pSync[SHMEM_REDUCE_SYNC_SIZE];
  static long pWrk[SHMEM_REDUCE_MIN_WRKDATA_SIZE];      

  for (int i = 0; i < SHMEM_REDUCE_SYNC_SIZE; i++)
    pSync[i] = SHMEM_SYNC_VALUE;
  for (int i = 0; i < SHMEM_REDUCE_MIN_WRKDATA_SIZE; i++)
    pWrk[i] = SHMEM_SYNC_VALUE;

  int npes = shmem_n_pes();
  int pe = shmem_my_pe();
  static long* PLOCKS = (long *) shmem_calloc(shmem_n_pes(), sizeof(long));                                     // Access to shared resources controlled by a single pe is determined by a lock on each pe
  static long* LOCK = (long *) shmem_calloc(1, sizeof(long)); 
  {
    Builder b(cli);
    Graph g = b.MakeGraph(pWrk, pSync);
    printf("Last check\n");
//    g.PrintTopology(LOCK);

/*    Partition<NodeID> vp(g.num_nodes());
    pvector<NodeID> comp(vp.max_width, true);
    if (pe == 0) {
      comp[0] = 0;
      comp[1] = 0;
      comp[2] = 4;
      comp[3] = 0;
      comp[4] = 4;
    } else {
      comp[0] = 4;
      comp[1] = 1;
      comp[2] = 5;
      comp[3] = 6;
      comp[4] = 5;
    }
    shmem_barrier_all();
    Compress(g, comp, vp);
    Link(6, 7, comp, vp);
    Compress(g, comp, vp);
    shmem_barrier_all();
    for (NodeID c : comp)
      printf("PE %d | c = %d\n", pe, c);
*/
/*    int64_t* thing = (int64_t*) shmem_calloc(5, sizeof(int64_t));
    for (int64_t i = 0; i < 5; i++)
      thing[i] = 4;
    shmem_barrier_all();
    if (pe == 0)
      if (shmem_int64_atomic_compare_swap(thing, 4, 15, 1))
        printf("Thank goodness, i was afraid c&s was broken\n");
*/
    auto CCBound = [](const Graph& gr){ return Afforest(gr, pSync); };
    //CCBound(g);
    BenchmarkKernel(cli, g, CCBound, PrintCompStats, CCVerifier);
  }

  shmem_finalize();
  return 0;
}
