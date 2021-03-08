#include <cinttypes>
#include <limits>
#include <iostream>
#include <queue>
#include <vector>

#include "benchmark.h"
#include "builder.h"
#include "command_line.h"
#include "graph.h"
#include "platform_atomics.h"
#include "pvector.h"
#include "timer.h"

using namespace std;

const long kDistInf = numeric_limits<WeightT>::max()/2;
const size_t kMaxBin = numeric_limits<size_t>::max()/2;
const size_t kBinSizeThreshold = 1000;

inline
void RelaxEdges(const WGraph &g, NodeID u, long delta, pvector<long> &dist, vector <vector<NodeID>> &local_bins, Partition<NodeID> vp) {
  long old_dist, new_dist;
  for (WNode wn : g.out_neigh(u)) {
    shmem_getmem(&old_dist, dist.begin()+vp.local_pos(wn.v), sizeof(long), vp.recv(wn.v));           
    shmem_getmem(&new_dist, dist.begin()+vp.local_pos(u), sizeof(long), vp.recv(u));               // can a pe fetch data from itself?
    new_dist += wn.w;
    while (new_dist < old_dist) {
      if (shmem_long_atomic_compare_swap(dist.begin()+vp.local_pos(wn.v), old_dist, new_dist, vp.recv(wn.v)) == old_dist) {              // not sure how to generically type atomics...
        size_t dest_bin = new_dist/delta;
        if (dest_bin >= local_bins.size())
          local_bins.resize(dest_bin+1);
        local_bins[dest_bin].push_back(wn.v);
        break;
      }
      shmem_getmem(&old_dist, dist.begin()+vp.local_pos(wn.v), sizeof(long), vp.recv(wn.v));               
    }
  }
}

pvector<long> Shmem_DeltaStep(const WGraph &g, NodeID source, int delta, long* pSync, long* pWrk) {
  long udist;
  long init = 0;
  long *old_dist, *new_dist;
  Timer t;
  printf("Delta: %d | source: %d\n", delta, source);
  Partition<NodeID> vp(g.num_nodes());                                  // Partition nodes: each PE is assigned a subset of nodes to process
  pvector<long> dist(vp.max_width, kDistInf, true);                     // Symmetric partitioned array of distances, initialized at infinity
  if (source >= vp.start && source < vp.end)                            // If src falls in your range of vertices, init distance to source as 0
    dist[vp.local_pos(source)] = init;
  Partition<NodeID> ep(g.num_edges_directed());                         // Partition edges: frontier is divided among PEs
  pvector<NodeID> frontier(g.num_edges_directed(), true);                 // Very, very bad! But the alternative is concurrent, resizable symmetric array 

  long* local_min = (long *) shmem_malloc(sizeof(long));
  size_t* frontier_tails = (size_t *) shmem_calloc(2, sizeof(size_t));        // init frontier tails to {1, 0}
  frontier_tails[0] = 1;
  long* shared_indexes = (long *) shmem_calloc(2, sizeof(long));        // size_t shared_indexes[2] = {0, kMaxBin};
  shared_indexes[1] = kMaxBin;
  frontier[0] = source;
  t.Start();                                                            // Timer start and stops are synch points
  vector<vector<NodeID> > local_bins(0);                                // Local vector of vector of node ids 
  NodeID* iter = (NodeID *) shmem_calloc(1, sizeof(NodeID));            // Shared counter for iteration (init 0)

  size_t copy_start;
  while (shared_indexes[(*iter)&1] != kMaxBin) {                        // if iter is odd, iter bitwise and (&) 1 is 1, if iter is even iter&1 is 0 and enter loop
    long &curr_bin_index = shared_indexes[(*iter)&1];                   // shared_indexes[0]    (always?)   
    long &next_bin_index = shared_indexes[((*iter)+1)&1];                      // shared_indexes[1]
    size_t &curr_frontier_tail = frontier_tails[(*iter)&1];                     // frontier_tails[0]
    size_t &next_frontier_tail = frontier_tails[((*iter)+1)&1];                 // frontier_tails[1]

    Partition<size_t> fp(curr_frontier_tail);
    for (size_t i = fp.start; i < fp.end; i++) {
      NodeID u = frontier[i];
      shmem_getmem(&udist, dist.begin()+vp.local_pos(u), sizeof(long), vp.recv(u));                      
      if (udist >= delta * static_cast<long>(curr_bin_index))
        RelaxEdges(g, u, delta, dist, local_bins, vp);
    }

    while (curr_bin_index < local_bins.size() &&
             !local_bins[curr_bin_index].empty() &&
             local_bins[curr_bin_index].size() < kBinSizeThreshold) {
      vector<NodeID> curr_bin_copy = local_bins[curr_bin_index];
      local_bins[curr_bin_index].resize(0);
      for (NodeID u : curr_bin_copy) {
        RelaxEdges(g, u, delta, dist, local_bins, vp);
      }
    }
    shmem_barrier_all();                                                // Is this needed? Next part of phase 1 is voting
    *local_min = next_bin_index;
    for (long i = curr_bin_index; i < local_bins.size(); i++) {         // Each PE finds its local min before voting
      if (!local_bins[i].empty()) {
        *local_min = min(*local_min, i);
        break;
      }
    }
    shmem_long_min_to_all(&next_bin_index, local_min, 1, 0, 0, vp.npes, pWrk, pSync); // next_bin_index = min of local_mins
    t.Stop();                                                          
    PrintStep(curr_bin_index, t.Millisecs(), curr_frontier_tail);       // End of phase 1

    t.Start();                                                                  // REMEMBER! Everyone must particpate in timer and label printing stuff
    curr_bin_index = kMaxBin;                                                   // Every PE updates current frontier tails and bin indexes to the same values
    curr_frontier_tail = 0;
    if (next_bin_index < local_bins.size()) {
      copy_start = shmem_ulong_atomic_fetch_add(&next_frontier_tail, local_bins[next_bin_index].size(), 0);       // PEs establish copying ranges using PE 0's copy of next frontier tail
      for (int i = 0; i < vp.npes; i++) {                                         // Copy your data into everyone's frontiers
        shmem_putmem(frontier.data() + copy_start, &*(local_bins[next_bin_index].begin()), sizeof(NodeID) * local_bins[next_bin_index].size(), i);
      }
      local_bins[next_bin_index].resize(0);
    }
    shmem_barrier_all();
    if (vp.pe == 0) {                                                           // couldn't every PE just update iter themselves?
      for (int i = 0; i < vp.npes; i++) {
        shmem_size_put(&next_frontier_tail, &next_frontier_tail, 1, i);
        shmem_int_atomic_inc(iter, i);
      }
    }
    shmem_barrier_all();
  }
  printf("PE %d | Took %d iterations\n", vp.pe, *iter);
  return dist;
}
 
void PrintSSSPStats(const WGraph &g, const pvector<long> &dist) {
  auto NotInf = [](long d) { return d != kDistInf; };
  int64_t num_reached = count_if(dist.begin(), dist.end(), NotInf);
  cout << "SSSP Tree reaches " << num_reached << " nodes" << endl;
}

// Print to file, compare results against original implementation
bool SSSPVerifier(const WGraph &g, NodeID source,
                  const pvector<long> &dist_to_test) {
  Partition<NodeID> vp(g.num_nodes());
  int* PRINTER = (int *) shmem_malloc(sizeof(int));
  *PRINTER = 0;
  shmem_barrier_all();
  shmem_int_wait_until(PRINTER, SHMEM_CMP_EQ, vp.pe);       // wait until previous PE puts your pe # in PRINTER
  ofstream shmem_out;
  shmem_out.open("/home/zach/projects/Dist_Mem_GAPBS/Dist_Mem_GAPBS/sssp_output.txt", ios::app);
  for (NodeID n = vp.start; n < vp.end; n++)
    shmem_out << dist_to_test[vp.local_pos(n)] << endl;
  shmem_out.close();
  if (!(vp.pe == vp.npes-1))
    shmem_int_p(PRINTER, vp.pe+1, vp.pe+1);             // who's next?
  /*if (vp.pe == 0) {
    ofstream shmem_out;
    shmem_out.open("/home/zach/projects/Dist_Mem_GAPBS/Dist_Mem_GAPBS/sssp_output.txt", ios::app);
    for (long dist : dist_to_test)
      shmem_out << dist << endl;
    shmem_out.close();
  }*/
  return true;
}


int main(int argc, char* argv[]) {
  CLDelta<WeightT> cli(argc, argv, "single-source shortest-path");
  if (!cli.ParseArgs())
    return -1;

  char size_env[] = "SMA_SYMMETRIC_SIZE=16G";
  putenv(size_env);

  shmem_init();

  static long pSync[SHMEM_REDUCE_SYNC_SIZE];
  static long pWrk[SHMEM_REDUCE_MIN_WRKDATA_SIZE];      
  static long long ll_pWrk[SHMEM_REDUCE_MIN_WRKDATA_SIZE];                 

  for (int i = 0; i < SHMEM_REDUCE_SYNC_SIZE; i++)
    pSync[i] = SHMEM_SYNC_VALUE;
  for (int i = 0; i < SHMEM_REDUCE_MIN_WRKDATA_SIZE; i++) {
    pWrk[i] = SHMEM_SYNC_VALUE;
    ll_pWrk[i] = SHMEM_SYNC_VALUE;
  }

  int pe = shmem_my_pe();
  int npes = shmem_n_pes();

  static long* PLOCKS = (long *) shmem_calloc(npes, sizeof(long));              // Access to shared resources controlled by a single pe is determined by a lock on each pe

  {
    WeightedBuilder b(cli);
    shmem_barrier_all();
    WGraph g = b.MakeGraph(pWrk, pSync);
    shmem_barrier_all();
    //g.PrintTopology();
    //g.PrintTopology(false);
    SourcePicker<WGraph> sp(g, cli.start_vertex());
    auto SSSPBound = [&sp, &cli] (const WGraph &g) {
      return Shmem_DeltaStep(g, sp.PickNext(), cli.delta(), pSync, pWrk);
    };
    SourcePicker<WGraph> vsp(g, cli.start_vertex());
    auto VerifierBound = [&vsp] (const WGraph &g, const pvector<long> &dist) {
      return SSSPVerifier(g, vsp.PickNext(), dist);
    };
    BenchmarkKernel(cli, g, SSSPBound, PrintSSSPStats, VerifierBound);
  }
  shmem_finalize();
  return 0;
}





