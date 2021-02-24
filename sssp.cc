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
      if (wn.v >= vp.start && wn.v < vp.end)
        old_dist = dist[vp.local_pos(wn.v)];      // swap failed, recheck dist update & retry
      else
        shmem_getmem(&old_dist, dist.begin()+vp.local_pos(wn.v), sizeof(long), vp.recv(wn.v));               
    }
  }
}


// Results: the symmetric distributions array is different on each PE
// distributions[i] = x means send x edges to PE i
// updates bin sizes and indices
void dist_scheme(long* distributions, long* bin_sizes, long* indices, Partition<NodeID> vp) {
  long m = 0;
  for (long i = 0; i < vp.npes; i++)
    m = max(m, bin_sizes[i]);
  long i = 0;
  while (i < vp.npes) {
    if (bin_sizes[i] == m)
      break;
    i++;
  }
  if (vp.pe == i) {                     // im the pe with the max number of elems
    for (long c = 0; c < vp.npes; c++) {
      if (bin_sizes[c] == 0 && bin_sizes[i] > 1) {
        distributions[c] = 1;
        //shmem_long_p(distributions+c, 1, c);
        bin_sizes[c]++;
        bin_sizes[i]--;
        shmem_long_atomic_inc(indices+c, c); 
        //indices[c]++;
      }
    }
  }
}

bool dist_heuristic(long* bin_sizes, int npes) {
  bool redistribute = false;
  long cumulative = 0;
  for (int i = 0; i < npes; i++) {                     // Sum local bin sizes
    cumulative += bin_sizes[i];
    if (bin_sizes[i] == 0) {
      redistribute = true;
    }
  }
  long avg = cumulative / npes;                      // some kind of heuristic to determine if work is balanced
  float dev = cumulative * (1.5 / npes);
  for (int i = 0; i < npes; i++) {
    if (bin_sizes[i] < (avg - dev) || bin_sizes[i] > (avg + dev)) {
      redistribute = true;
      break;
    }
  }
  if (cumulative < npes)
    redistribute = false;
 /* if (redistribute) {
    for (int i = 0; i < npes; i++)
      shmem_int_p(redist, 1, i);
  }*/
  printf("Cumulative: %lu | redist: %d\n", cumulative, (int) redistribute);
  return(redistribute);
}

// If each PE has a reasonable bin size, just copy local bins into local frontier
// Otherwise distribute bin contents among PEs
// Calculate the new frontier size for each PE and return it
long distribute_frontier(pvector<long> &frontier, vector <vector<NodeID>> &local_bins, 
                  long* bin_sizes, long next_bin_index, long next_frontier_tail, Partition<NodeID> vp, Partition<NodeID> ep, long* indices, long* distributions) {
  long return_val;
  static long local_pSync[SHMEM_BCAST_SYNC_SIZE];
  for (int i = 0; i < SHMEM_BCAST_SYNC_SIZE; i++)
    local_pSync[i] = SHMEM_SYNC_VALUE;

  //long* indices = (long *) shmem_calloc(vp.npes, sizeof(long)); // Number of elements you recieved from other PEs in dist_scheme
  // although really it makes more sense for indices[i] = x on PE k to mean k recieves x edges from PE i (notes for advanced distribution scheme)
  //long* distributions = (long *) shmem_calloc(vp.npes, sizeof(long)); // Each PE maintains a list: distributions[i] = x means send x edges to PE i

  // Reset the indices/distributions arrays
  for (int i = 0; i < vp.npes; i++) {
    indices[i] = 0;
    distributions[i] = 0;
  }
  shmem_barrier_all();

  // Determine if work must be re-balanced
  int* redistribute = (int*) shmem_calloc(1, sizeof(int));
  if (vp.pe == 0) {
    if (dist_heuristic(bin_sizes, vp.npes)) { 
      for (int i = 0; i < vp.npes; i++)
        shmem_int_p(redistribute, 1, i);
    }
  }
  shmem_barrier_all();
  //printf("PE %d - Redistribute: %d\n", vp.pe, *redistribute);

  // Update all PEs with PE 0's copy of bin_sizes
  shmem_broadcast64(bin_sizes, bin_sizes, vp.npes, 0, 0, 0, vp.npes, local_pSync);
  shmem_barrier_all();

  if (*redistribute == 1) {                                   // Need a better distribution scheme, right now i just want everyone to have an edge
    //printf("redisting\n");
    dist_scheme(distributions, bin_sizes, indices, vp);
    shmem_barrier_all();
    long count = 0;
    for (long i = 0; i < vp.npes; i++) {
      if (distributions[i] > 0) {
        //printf("PE %d wants to put an edge to PE %d\n", vp.pe, i);
        shmem_long_p(frontier.data(), *(local_bins[next_bin_index].begin()+count), i);
        count++;
      }
    }
    //printf("2.3\n");
    shmem_barrier_all();
    long binsize = 0;
    if (next_bin_index < local_bins.size()) {
      binsize = local_bins[next_bin_index].size();
      //printf("check 2.4\n");
      copy(local_bins[next_bin_index].begin()+count, local_bins[next_bin_index].end(), frontier.data()+indices[vp.pe]); // copy the edges you didn't give away into frontier starting at how many edges you received
      local_bins[next_bin_index].resize(0);
    }
    //printf("check 2.3.1\n");
    return_val = indices[vp.pe] + (binsize-count); 
    /*shmem_barrier_all();
    shmem_free(indices);
    shmem_free(distributions);
    shmem_free(redistribute);*/
    return return_val;
  } else {                                              // Copy contents of local bin into personal frontier
    //printf("check 4\n");
    if (next_bin_index < local_bins.size()) {
      long binsize = local_bins[next_bin_index].size();
      copy(local_bins[next_bin_index].begin(), local_bins[next_bin_index].end(), frontier.data());
      local_bins[next_bin_index].resize(0);
      return(binsize);
    } else {
      return(next_frontier_tail);
      //return(next_bin_index);
      //return 0;
    }
  }
} 

pvector<long> Shmem_DeltaStep(const WGraph &g, NodeID source, int delta, long* VOTE_LOCK, long* SINGLE_LOCK, long* pSync, long* pWrk) {
  long udist;
  long init = 0;
  long *old_dist, *new_dist;
  Timer t;
  Partition<NodeID> vp(g.num_nodes());                                  // Partition nodes: each PE is assigned a subset of nodes to process
  pvector<long> dist(vp.max_width, kDistInf, true);                     // Symmetric partitioned array of distances, initialized at infinity
  if (source >= vp.start && source < vp.end)                            // If src falls in your range of vertices, init distance to source as 0
    dist[vp.local_pos(source)] = init;
  Partition<NodeID> ep(g.num_edges_directed());                         // Partition edges: frontier is divided among PEs
  pvector<long> frontier(ep.max_width, true);                         // Symmetric partitioned frontier

  long* indices = (long *) shmem_calloc(vp.npes, sizeof(long)); // Number of elements you recieved from other PEs in dist_scheme
  long* distributions = (long *) shmem_calloc(vp.npes, sizeof(long)); // Each PE maintains a list: distributions[i] = x means send x edges to PE i

  long* local_min = (long *) shmem_malloc(sizeof(long));                // Allocate voting space for next bin at the end of phase 1
  long* bin_sizes = (long *) shmem_calloc(vp.npes, sizeof(long));        // PE 0 maintains a bin size for each PE to update at the end of each iteration
  long* frontier_tails = (long *) shmem_calloc(2, sizeof(long));        // init frontier tails to {1, 0}
  frontier_tails[0] = 1;
  long* shared_indexes = (long *) shmem_calloc(2, sizeof(long));        // size_t shared_indexes[2] = {0, kMaxBin};
  shared_indexes[1] = kMaxBin;
  if (vp.pe == 0) {                                                     // PE 0 contains src as the first frontier element
    frontier[0] = source;
  }
  t.Start();                                                            // Timer start and stops are synch points
  vector<vector<NodeID> > local_bins(0);                                // Local vector of vector of node ids 
  NodeID* iter = (NodeID *) shmem_calloc(1, sizeof(NodeID));            // Shared counter for iteration (init 0)

  while (shared_indexes[(*iter)&1] != kMaxBin) {                        // if iter is odd, iter bitwise and (&) 1 is 1, if iter is even iter&1 is 0 and enter loop
    long &curr_bin_index = shared_indexes[(*iter)&1];                   // shared_indexes[0]    (always?)   
    long &next_bin_index = shared_indexes[((*iter)+1)&1];                      // shared_indexes[1]
    long &curr_frontier_tail = frontier_tails[(*iter)&1];                     // frontier_tails[0]
    long &next_frontier_tail = frontier_tails[((*iter)+1)&1];                 // frontier_tails[1]

    for (long i = 0; i < curr_frontier_tail; i++) {                     // Each PE processes their own frontier (asynchronous)
      NodeID u = frontier[i];
      if (u >= vp.start && u < vp.end)
        udist = dist[vp.local_pos(u)];
      else
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
    for (long i=curr_bin_index; i < local_bins.size(); i++) {         // Each PE finds its local min before voting
      if (!local_bins[i].empty()) {
        *local_min = min(*local_min, i);
        break;
      }
    }
    //shmem_long_min_to_all(shared_indexes+(((*iter)+1)&1), local_min, 1, 0, 0, vp.npes, pWrk, pSync); // next_bin_index = min of local_mins
    shmem_long_min_to_all(&next_bin_index, local_min, 1, 0, 0, vp.npes, pWrk, pSync); // next_bin_index = min of local_mins
    t.Stop();                                                          
    PrintStep(curr_bin_index, t.Millisecs(), curr_frontier_tail);       // End of phase 1

    t.Start();                                                                  // REMEMBER! Everyone must particpate in timer and label printing stuff
    curr_bin_index = kMaxBin;                                                   // Every PE updates current frontier tails and bin indexes to the same values
    curr_frontier_tail = 0;
    shmem_barrier_all();

    if (next_bin_index < local_bins.size()) {
      shmem_long_p(bin_sizes+vp.pe, local_bins[next_bin_index].size(), 0);      // Add your local bin size to the stable array on PE 0
    } else {
      shmem_long_p(bin_sizes+vp.pe, 0, 0);
    }
    shmem_barrier_all();
    next_frontier_tail = distribute_frontier(frontier, local_bins, bin_sizes, next_bin_index, next_frontier_tail, vp, ep, indices, distributions);       // Each PE gets their own frontier tail representing size of local frontier
    //(*iter)++;
    shmem_barrier_all();
    if (vp.pe == 0) {                                                           // couldn't every PE just update iter themselves?
      for (int i = 0; i < vp.npes; i++) {
        shmem_int_atomic_inc(iter, i);
      }
    }
    //printf("PE %d | nft = %lu\n", vp.pe, next_frontier_tail);
    //for (long l = 0; l < next_frontier_tail; l++)
      //printf("PE %d | Frontier[%lu] = %lu\n", vp.pe, l, frontier[l]);
    shmem_barrier_all();
  }
  printf("PE %d | Took %d iterations\n", vp.pe, *iter);
  /*shmem_barrier_all();
  pvector<long> combined_dists = dist.combine(g.num_nodes(), pSync);
  return combined_dists;*/
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
    auto discard = sp.PickNext();
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





