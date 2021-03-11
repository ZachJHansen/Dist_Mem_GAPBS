// Copyright (c) 2015, The Regents of the University of California (Regents)
// See LICENSE.txt for license details

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


/*
GAP Benchmark Suite
Kernel: Single-source Shortest Paths (SSSP)
Author: Scott Beamer

Returns array of distances for all vertices from given source vertex

This SSSP implementation makes use of the ∆-stepping algorithm [1]. The type
used for weights and distances (WeightT) is typedefined in benchmark.h. The
delta parameter (-d) should be set for each input graph.

The bins of width delta are actually all thread-local and of type std::vector
so they can grow but are otherwise capacity-proportional. Each iteration is
done in two phases separated by barriers. In the first phase, the current
shared bin is processed by all threads. As they find vertices whose distance
they are able to improve, they add them to their thread-local bins. During this
phase, each thread also votes on what the next bin should be (smallest
non-empty bin). In the next phase, each thread copies their selected
thread-local bin into the shared bin.

Once a vertex is added to a bin, it is not removed, even if its distance is
later updated and it now appears in a lower bin. We find ignoring vertices if
their current distance is less than the min distance for the bin to remove
enough redundant work that this is faster than removing the vertex from older
bins.

[1] Ulrich Meyer and Peter Sanders. "δ-stepping: a parallelizable shortest path
    algorithm." Journal of Algorithms, 49(1):114–152, 2003.
*/



using namespace std;

const long kDistInf = numeric_limits<WeightT>::max()/2;
const size_t kMaxBin = numeric_limits<size_t>::max()/2;
const size_t kBinSizeThreshold = 1000;

inline
void RelaxEdges(const WGraph &g, NodeID u, long delta, pvector<long> &dist, vector <vector<NodeID>> &local_bins, Partition<NodeID> vp) {
  long old_dist, new_dist;
  //printf("PE %d relaxing edge %d\n", vp.pe, u);
  for (WNode wn : g.out_neigh(u)) {
    // modifications to local data should also be atomic! 
    shmem_getmem(&old_dist, dist.begin()+vp.local_pos(wn.v), sizeof(long), vp.recv(wn.v));           
    shmem_getmem(&new_dist, dist.begin()+vp.local_pos(u), sizeof(long), vp.recv(u));               // can a pe fetch data from itself?
    new_dist += wn.w;
    while (new_dist < old_dist) {
      // orig imp returns true if old_dist was the same as the elem at dist[wn.v]
      // but shmem C&S returns the value at dist[wn.v], not a bool
      //if (shmem_long_atomic_compare_swap(dist.begin()+vp.local_pos(wn.v), old_dist, new_dist, vp.recv(wn.v))) {              // not sure how to generically type atomics...
      if (shmem_long_atomic_compare_swap(dist.begin()+vp.local_pos(wn.v), old_dist, new_dist, vp.recv(wn.v)) == old_dist) {              // not sure how to generically type atomics...
        size_t dest_bin = new_dist/delta;
        if (dest_bin >= local_bins.size())
          local_bins.resize(dest_bin+1);
        local_bins[dest_bin].push_back(wn.v);
        break;
      }
      if (wn.v >= vp.start && wn.v < vp.end)
        old_dist = dist[wn.v];      // swap failed, recheck dist update & retry
      else
        shmem_getmem(&old_dist, dist.begin()+vp.local_pos(wn.v), sizeof(long), vp.recv(wn.v));               
    }
  }
}

// should return pvector<WeightT> but thats a v2 thing
pvector<long> Shmem_DeltaStep(const WGraph &g, NodeID source, int delta, long* VOTE_LOCK, long* SINGLE_LOCK, long* pSync, long* pWrk) {
  long udist;
  long init = 0;
  long *old_dist, *new_dist;
  Timer t;
  //printf("Source: %d\n", source);
  Partition<NodeID> vp(g.num_nodes());
//  pvector<WeightT> dist(vp.max_width, kDistInf, true);                                 // Initialize symmetric partitioned pvector with +inf
  pvector<long> dist(vp.max_width, kDistInf, true);
  if (source >= vp.start && source < vp.end)                            // If src falls in your range of vertices, init distance to source as 0
    dist[vp.local_pos(source)] = init;
  //else
    //shmem_putmem(dist.begin()+vp.local_pos(source), &init, sizeof(long), vp.recv(source));
  //printf("PE %d is init frontier with size = %lu\n", vp.pe, g.num_edges_directed());

  // we could probably partition frontier, would add communication overhead
  pvector<NodeID> frontier(g.num_edges_directed(), true);                               // Symmetric frontier - intentionally overallocated. Do we want a space or performance hit?
  // two element arrays for double buffering curr=iter&1, next=(iter+1)&1
  long* shared_indexes = (long *) shmem_calloc(2, sizeof(long));                  // size_t shared_indexes[2] = {0, kMaxBin};
  shared_indexes[1] = kMaxBin;
  long* frontier_tails = (long *) shmem_calloc(2, sizeof(long));                   // size_t frontier_tails[2] = {1, 0};
  frontier_tails[0] = 1;
  frontier[0] = source;
  t.Start();                                                                    // Timer start and stops are synch points
  vector<vector<NodeID> > local_bins(0);                                      // vector of vector of node ids (thread local?)
  NodeID* iter = (NodeID *) shmem_calloc(1, sizeof(NodeID));                                  // Shared counter for next node (init 0)
  //printf("PE %d | dist.begin = %p, dist.end = %p, front.begin = %p, front.end = %p\n", vp.pe, (void *) dist.begin(),  (void *) dist.end(), (void *) frontier.begin(), (void *) frontier.end());
  while (shared_indexes[(*iter)&1] != kMaxBin) {                                 // if iter is odd, iter bitwise and (&) 1 is 1, if iter is even iter&1 is 0 and enter loop
    //printf("PE %d | While start | iter = %d, curr_front = %lu\n", vp.pe, *iter, frontier_tails[(*iter)&1]);
    // these all reference the shared data structs
    long &curr_bin_index = shared_indexes[(*iter)&1];                          
    long &next_bin_index = shared_indexes[((*iter)+1)&1];                      
    long &curr_frontier_tail = frontier_tails[(*iter)&1];
    long &next_frontier_tail = frontier_tails[((*iter)+1)&1];
    //*aye = 0;
    //printf("PE %d | post init | iter = %d, curr_front = %lu\n", vp.pe, *iter, curr_frontier_tail);
    //printf("Pe %d | While start | iter = %d, shared_index[iter&1] = %lu\n", vp.pe, *iter, shared_indexes[(*iter)&1]);
    Partition<NodeID> fp(curr_frontier_tail);                                            // I don't know how to do a dynamic schedule nowait, so I'm just gonna naively partition the work
    //printf("Everyone here? PE %d\n", vp.pe);
    //shmem_barrier_all();
    //printf("PE %d fp.start: %lu | fp.end: %lu\n", fp.pe, fp.start, fp.end);
    //printf("PE %d | Current frontier tail: %lu\n", vp.pe, frontier_tails[(*iter)&1]);


    // maybe, instead of naively partitioning the work, have each pe loop through the entire frontier
    // if a PE finds a node assigned to them they process it otherwise leave it alone
    for (long i = fp.start; i < fp.end; i++) {
      NodeID u = frontier[i];
      if (u >= vp.start && u < vp.end)
        //udist = dist[u];
        udist = dist[vp.local_pos(u)];
      else
        shmem_getmem(&udist, dist.begin()+vp.local_pos(u), sizeof(long), vp.recv(u));                       
      //printf("i: %d | Udist: %d\n", i, udist);
      //printf("dist[%d] = %d, curr_bin_index = %d, delta = %d\n", u, udist, static_cast<WeightT>(curr_bin_index), delta);
      if (udist >= delta * static_cast<long>(curr_bin_index))
        //printf("the check\n");
        RelaxEdges(g, u, delta, dist, local_bins, vp);
    }
    while (curr_bin_index < local_bins.size() &&
             !local_bins[curr_bin_index].empty() &&
             local_bins[curr_bin_index].size() < kBinSizeThreshold) {
      //printf("mystery\n");
      vector<NodeID> curr_bin_copy = local_bins[curr_bin_index];
      local_bins[curr_bin_index].resize(0);
      for (NodeID u : curr_bin_copy) {
        //printf("Relaxing edge: %d\n", u);
        RelaxEdges(g, u, delta, dist, local_bins, vp);
      }
    }
    shmem_barrier_all();
    //shmem_global_exit(0);
    //exit(0);

    long* local_min = (long *) shmem_malloc(sizeof(long));        // maybe alloc this before while loop 
    *local_min = next_bin_index;
    for (long i=curr_bin_index; i < local_bins.size(); i++) {         // Each PE finds its local min before voting
      if (!local_bins[i].empty()) {
        *local_min = min(*local_min, i);
        break;
      }
    }
    //shmem_long_min_to_all(shared_indexes+(((*iter)+1)&1), local_min, 1, 0, 0, vp.npes, pWrk, pSync); // next_bin_index = min of local_mins
    shmem_long_min_to_all(&next_bin_index, local_min, 1, 0, 0, vp.npes, pWrk, pSync); // next_bin_index = min of local_mins
    //if (curr_bin_index < local_bins.size())
      //printf("PE %d check 1\n", vp.pe);
    /*for (size_t i=curr_bin_index; i < local_bins.size(); i++) {         // this region could def be optimized, maybe each pe finds its local min b4 voting?
      if (!local_bins[i].empty()) {
        shmem_set_lock(VOTE_LOCK);                                            // critical region
        for (int p = 0; p < vp.npes; p++) {
          size_t temp = min(shared_indexes[((*iter)+1)&1], i);
          //printf("PE %d is voting: temp = %lu\n", vp.pe, temp);
          shmem_putmem(shared_indexes+(((*iter)+1)&1), &temp, sizeof(size_t), p);
        }
        shmem_clear_lock(VOTE_LOCK);
        break;
      }      
    }*/
    shmem_barrier_all();
    printf("PE %d | next chosen bin is %lu\n", vp.pe, shared_indexes[((*iter)+1)&1]);
    //shmem_global_exit(0);
    //exit(0);
    t.Stop();
    PrintStep(curr_bin_index, t.Millisecs(), curr_frontier_tail);
    t.Start();                                                          //REMEMBER! Everyone must particpate in timer and label printing stuff
    //if (shmem_test_lock(SINGLE_LOCK) == 0) {                                  // lock is not set, execute body 
    curr_bin_index = kMaxBin;                                            // so since &curr_bin_index = shared_index, do changes to curr_bin_index affect (or should affect) the shared?
    //shared_indexes[(*iter)&1] = kMaxBin;
    curr_frontier_tail = 0;
    //frontier_tails[(*iter)&1] = 0;
    //printf("PE %d has set curr bin to kmax and curr front to 0\n", vp.pe);
    //printf("PE %d | Curr bin index = %lu | curr_front = %lu | shared_indexes[iter&1] = %lu | shared_indexes[(iter+1)&1] = %lu\n", vp.pe, curr_bin_index, curr_frontier_tail, shared_indexes[(*iter)&1], shared_indexes[((*iter)+1)&1]);


    //shmem_clear_lock(SINGLE_LOCK);
    //}
    //printf("PE %d checking in\n", vp.pe);
    //printf("PE %d | iter: %d \n", vp.pe, *iter);
    if (next_bin_index < local_bins.size()) {
      printf("PE %d Frontier size = %lu\n", vp.pe, frontier.size());
      printf("Frontier contents b4 copy: ");
      for (WeightT t : frontier)
        printf("%d ", t);
      printf("\n");
      shmem_set_lock(SINGLE_LOCK);
      size_t copy_start;
      for (int i = 0; i < vp.npes; i++) {
        //if (vp.pe != i) {
        //copy_start = shmem_ulong_atomic_fetch_add(&next_frontier_tail, local_bins[next_bin_index].size(), i);
        copy_start = shmem_long_atomic_fetch_add(&next_frontier_tail, local_bins[next_bin_index].size(), i);
//        printf("Copy start: %lu\n", copy_start);
        //printf("local.begin = %d\n", *(local_bins[next_bin_index].begin()));
 //       printf("end - begin = %d | size = %d\n", local_bins[next_bin_index].end() - local_bins[next_bin_index].begin(), local_bins[next_bin_index].size());
        shmem_putmem(frontier.data() + copy_start, (void*) local_bins[next_bin_index].data(), sizeof(NodeID) * local_bins[next_bin_index].size(), i);
      }
  /*    printf("PE %d | Frontier contents after copy: ", vp.pe);
      for (WeightT t : frontier)
        printf("%d ", t);
      printf("\n");
*/
      shmem_clear_lock(SINGLE_LOCK);
      local_bins[next_bin_index].resize(0);
    }

    //printf("Pe %d | Prior to iter inc | iter = %d, shared_index[iter&1] = %lu, curr_front = %lu\n", vp.pe, *iter, shared_indexes[(*iter)&1], curr_frontier_tail);
    //iter++;
    if (vp.pe == 0) {
      for (int i = 0; i < vp.npes; i++) {
        shmem_int_atomic_inc(iter, i);
        //printf("PE %d just incremented iter on PE %d\n", vp.pe, i);
      }
    }
    shmem_barrier_all();
    //printf("PE %d | While end | curr_front = %lu, next_front = %lu\n", vp.pe, 
    //printf("Pe %d | While end | iter = %d, shared_index[iter&1] = %lu, curr_front = %lu\n", vp.pe, *iter, shared_indexes[(*iter)&1], curr_frontier_tail);
    //printf("PE %d | While end | iter = %d, curr_frontier_tail = %lu\n", vp.pe, *iter, frontier_tails[(*iter)&1]); //curr_frontier_tail);
    //printf("Pe %d | After iter inc | iter = %d, shared_index[iter&1] = %lu, curr_front = %lu\n", vp.pe, *iter, shared_indexes[(*iter)&1], curr_frontier_tail);
  }

  printf("PE %d | Took %d iterations\n", vp.pe, *iter);
  //printf("PE %d is here\n", vp.pe);
  /*if (shmem_test_lock(SINGLE_LOCK) == 0) {
    //std::cout << "took " << iter << " iterations" << endl;
    printf("PE %d | Took %d iterations\n", vp.pe, *iter);
    shmem_clear_lock(SINGLE_LOCK);
  }*/
  shmem_barrier_all();
  pvector<long> combined_dists = dist.combine(g.num_nodes(), pSync);
  //shmem_collect32(combined_dists.begin(), dist.begin(), vp.end-vp.start, 0, 0, vp.npes, pSync);
  return combined_dists;
  //return dist;                                                          // only pe 0 has a correct distance array
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
  if (vp.pe == 0) {
    ofstream shmem_out;
    shmem_out.open("/home/zach/projects/Dist_Mem_GAPBS/Dist_Mem_GAPBS/sssp_output.txt", ios::app);
    for (long dist : dist_to_test)
      shmem_out << dist << endl;
    shmem_out.close();
  }
}


int main(int argc, char* argv[]) {
  CLDelta<WeightT> cli(argc, argv, "single-source shortest-path");
  if (!cli.ParseArgs())
    return -1;

  char size_env[] = "SMA_SYMMETRIC_SIZE=16G";
  putenv(size_env);

  static long PRINT_LOCK = 0;
  static long VOTE_LOCK = 0;
  static long SINGLE_LOCK = 0;

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
    //printf("Check 1\n");
    shmem_barrier_all();
    WGraph g = b.MakeGraph(pWrk, pSync);
    //printf("Check 2\n");
    shmem_barrier_all();
    //g.PrintTopology();
    //g.PrintTopology(false);
    SourcePicker<WGraph> sp(g, cli.start_vertex());
    /*Partition vp(g.num_nodes());
    printf("PE %d is init dist with width = %lu\n", vp.pe, vp.max_width);
    pvector<WeightT> dist(vp.max_width, kDistInf, true);                                 // Initialize symmetric partitioned pvector with +inf
    WeightT init = 0;
    NodeID source = 1;
    if (source >= vp.start && source < vp.end)
      dist[vp.local_pos(source)] = init;
    else
      shmem_putmem(dist.begin()+vp.local_pos(source), &init, sizeof(WeightT), vp.recv(source));
    printf("PE %d is init frontier with size = %lu\n", vp.pe, g.num_edges_directed());
    pvector<NodeID> frontier(g.num_edges_directed(), true);  
    if (shmem_int_atomic_compare_swap(dist.begin(), kDistInf, 15, 0))
      printf("swapped\n");*/
    //pvector<WeightT> parent = SSSPBound(g);*/
    /*shmem_barrier_all();
    if (shmem_my_pe() == 0) {
      int i = 0;        
      for (long p : parent)
        printf("PE %d  | dist(%d) = %lu\n", shmem_my_pe(), i++, p);
    }*/
    //printf("escape\n");  
    auto SSSPBound = [&sp, &cli] (const WGraph &g) {
      return Shmem_DeltaStep(g, sp.PickNext(), cli.delta(), &VOTE_LOCK, &SINGLE_LOCK, pSync, pWrk);
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





