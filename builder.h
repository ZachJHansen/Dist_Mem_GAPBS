// Copyright (c) 2015, The Regents of the University of California (Regents)
// See LICENSE.txt for license details

#ifndef BUILDER_H_
#define BUILDER_H_

#include <stdlib.h>
#include <algorithm>
#include <cinttypes>
#include <fstream>
#include <functional>
#include <type_traits>
#include <utility>

#include "tournament.h"
#include "command_line.h"
#include "generator.h"
#include "graph.h"
#include "platform_atomics.h"
#include "pvector.h"
#include "reader.h"
#include "timer.h"
#include "util.h"


/*
GAP Benchmark Suite
Class:  BuilderBase
Author: Scott Beamer
Given arguements from the command line (cli), returns a built graph
 - MakeGraph() will parse cli and obtain edgelist and call
   MakeGraphFromEL(edgelist) to perform actual graph construction
 - edgelist can be from file (reader) or synthetically generated (generator)
 - Common case: BuilderBase typedef'd (w/ params) to be Builder (benchmark.h)
*/


template <typename NodeID_, typename DestID_ = NodeID_,
          typename WeightT_ = NodeID_, bool invert = true>
class BuilderBase {
  typedef EdgePair<NodeID_, DestID_> Edge;
  typedef pvector<Edge> EdgeList;

  const CLBase &cli_;
  bool symmetrize_;
  bool needs_weights_;
  int64_t num_nodes_ = -1;

 public:
  explicit BuilderBase(const CLBase &cli) : cli_(cli) {
    symmetrize_ = cli_.symmetrize();
    needs_weights_ = !std::is_same<NodeID_, DestID_>::value;
  }

  DestID_ GetSource(EdgePair<NodeID_, NodeID_> e) {
    return e.u;
  }

  DestID_ GetSource(EdgePair<NodeID_, NodeWeight<NodeID_, WeightT_>> e) {
    return NodeWeight<NodeID_, WeightT_>(e.u, e.v.w);
  }

  NodeID_ FindMaxNodeID(const EdgeList &el) {
    long* pSync = (long*) shmem_calloc(SHMEM_REDUCE_SYNC_SIZE, sizeof(long));
    int* pWrk = (int*) shmem_calloc(SHMEM_REDUCE_MIN_WRKDATA_SIZE, sizeof(int));
    for (int i = 0; i < SHMEM_REDUCE_SYNC_SIZE; i++)
      pSync[i] = SHMEM_SYNC_VALUE;
    for (int i = 0; i < SHMEM_REDUCE_MIN_WRKDATA_SIZE; i++) 
      pWrk[i] = SHMEM_SYNC_VALUE;
    shmem_barrier_all();                                      // not needed due to calloc?
    NodeID_* max_seen = (NodeID_*) shmem_calloc(1, sizeof(NodeID_));
    for (auto it = el.begin(); it < el.end(); it++) {              // find local max  - is el.end() max width (possibly uninit values) or local width?
      Edge e = *it;
      *max_seen = std::max(*max_seen, e.u);
      *max_seen = std::max(*max_seen, (NodeID_) e.v);
    }
    shmem_barrier_all();
    shmem_int_max_to_all(max_seen, max_seen, 1, 0, 0, shmem_n_pes(), pWrk, pSync);
    return *max_seen;
  }

  // Return pvector representing degrees for vertices assigned to local PE
  // pvector is symmetric and up-to-date but unsynched - do not synch!
  // pvectors on each pe should be concatenated to make a complete list (once unused remainder is trimmed off PEs != npes-1)
  // ASSUMES NODEIDS ARE INTS
  pvector<NodeID_> CountDegrees(const EdgeList &el, bool transpose, Partition<NodeID_>* vp, Partition<>* ep = NULL) {
    int local_v, receiver;
    pvector<NodeID_> degrees(vp->max_width, 0, true);                                     // Symmetric pvector of size max partition width
    //#pragma omp parallel for
    Edge e;
    int flush_counter = 0;
    shmem_barrier_all();
    for (auto it = el.begin(); it < el.end(); it++) {
      e = *it;
      if (symmetrize_ || (!symmetrize_ && !transpose)) {
        receiver = vp->recv(e.u);
        local_v = vp->local_pos(e.u);
        shmem_int_atomic_inc(degrees.begin()+(local_v), receiver);                                   // increment degree of vertex e.u on pe receiver (could be local PE) 
      }
      if (symmetrize_ || (!symmetrize_ && transpose)) {
        receiver = vp->recv(e.v);
        local_v = vp->local_pos(e.v);
        shmem_int_atomic_inc(degrees.begin()+(local_v), receiver);                                   // increment degree of vertex e.v on pe receiver (could be local PE) 
      }
      flush_counter++;
      if (flush_counter % 2000000 == 0)  // weirdness: without periodic barriers, CountDegrees runs out of memory on twitter, road. barrier forces shmem to flush communication buffers maybe?
        shmem_barrier_all();
    }
    shmem_barrier_all();
    return degrees;
  }

  static
  pvector<SGOffset> PrefixSum(const pvector<NodeID_> &degrees) {
    pvector<SGOffset> sums(degrees.size() + 1);
    SGOffset total = 0;
    for (size_t n=0; n < degrees.size(); n++) {
      sums[n] = total;
      total += degrees[n];
    }
    sums[degrees.size()] = total;
    return sums;
  }

  static
  pvector<SGOffset> ParallelPrefixSum(const pvector<NodeID_> &degrees) {
    const size_t block_size = 1<<20;
    const size_t num_blocks = (degrees.size() + block_size - 1) / block_size;
    pvector<SGOffset> local_sums(num_blocks);
    #pragma omp parallel for
    for (size_t block=0; block < num_blocks; block++) {
      SGOffset lsum = 0;
      size_t block_end = std::min((block + 1) * block_size, degrees.size());
      for (size_t i=block * block_size; i < block_end; i++)
        lsum += degrees[i];
      local_sums[block] = lsum;
    }
    pvector<SGOffset> bulk_prefix(num_blocks+1);
    SGOffset total = 0;
    for (size_t block=0; block < num_blocks; block++) {
      bulk_prefix[block] = total;
      total += local_sums[block];
    }
    bulk_prefix[num_blocks] = total;
    pvector<SGOffset> prefix(degrees.size() + 1, true);
    #pragma omp parallel for
    for (size_t block=0; block < num_blocks; block++) {
      SGOffset local_total = bulk_prefix[block];
      size_t block_end = std::min((block + 1) * block_size, degrees.size());
      for (size_t i=block * block_size; i < block_end; i++) {
        prefix[i] = local_total;
        local_total += degrees[i];
      }
    }
    prefix[degrees.size()] = bulk_prefix[num_blocks];
    return prefix;
  }

  // Removes self-loops and redundant edges
  // Side effect: neighbor IDs will be sorted
  void SquishCSR(const CSRGraph<NodeID_, DestID_, invert> &g, bool transpose,
                 DestID_*** sq_index, DestID_** sq_neighs, Partition<NodeID_> vp, long* pSync, long* pWrk) {
    NodeID_ indx;
    pvector<NodeID_> diffs(vp.max_width);
    DestID_ *n_start, *n_end;
    #pragma omp parallel for private(n_start, n_end)
    for (NodeID_ n = vp.start; n < vp.end; n++) {
      indx = vp.local_pos(n);
      if (transpose) {
        n_start = g.in_neigh(n).start();
        n_end = g.in_neigh(n).finish();
        //printf("PE %d | n: %d | n_start: %p => %d n_end: %p => %d (t)\n", vp.pe, n, (void*) n_start, *n_start, (void*) n_end, *n_end);
      } else {
        n_start = g.out_neigh(n).start();
        n_end = g.out_neigh(n).finish();
        //printf("PE %d | n: %d | n_start: %p => %d n_end: %p => %d (nt)\n", vp.pe, n, (void*) n_start, *n_start, (void*) n_end, *n_end);
      }
      std::sort(n_start, n_end);                        // sorts elements in range [n_start, n_end)
      DestID_ *new_end = std::unique(n_start, n_end);   // removes all but one of consecutive elements
      new_end = std::remove(n_start, new_end, n);       // removes self-loops (v,v)
      diffs[indx] = new_end - n_start;
    }
    pvector<SGOffset> sq_offsets = ParallelPrefixSum(diffs);
    SGOffset* max_neigh = (SGOffset *) shmem_malloc(sizeof(SGOffset));
    shmem_long_max_to_all(max_neigh, sq_offsets.begin()+(vp.end - vp.start), 1, 0, 0, vp.npes, pWrk, pSync); 
    *sq_neighs = (DestID_ *) shmem_calloc(*max_neigh, sizeof(DestID_));
    *sq_index = CSRGraph<NodeID_, DestID_>::GenIndex(sq_offsets, *sq_neighs, &vp);
    shmem_barrier_all();
    #pragma omp parallel for private(n_start)
    for (NodeID_ n=vp.start; n < vp.end; n++) {
      indx = vp.local_pos(n);
      if (transpose)
        n_start = g.in_neigh(n).start();        // transpose for incoming neighbors
      else
        n_start = g.out_neigh(n).start();
      std::copy(n_start, n_start+diffs[indx], (*sq_index)[indx]);    // copy elements [n_start, n_start+diffs[indx]) into (*sq_index)[indx]
    }
    shmem_barrier_all();
  }

  CSRGraph<NodeID_, DestID_, invert> SquishGraph(
      const CSRGraph<NodeID_, DestID_, invert> &g, Partition<NodeID_>* vp, long* pSync, long* pWrk) {
    DestID_ **out_index, *out_neighs, **in_index, *in_neighs;
    SquishCSR(g, false, &out_index, &out_neighs, *vp, pSync, pWrk);
    shmem_barrier_all();
    if (g.directed()) {
      shmem_barrier_all();
      if (invert)
        SquishCSR(g, true, &in_index, &in_neighs, *vp, pSync, pWrk);
      return CSRGraph<NodeID_, DestID_, invert>(g.num_nodes(), out_index,
                                                out_neighs, in_index,
                                                in_neighs, pSync, pWrk);
    } else {
      shmem_barrier_all();
      return CSRGraph<NodeID_, DestID_, invert>(g.num_nodes(), out_index,
                                                out_neighs, pSync, pWrk);
    }
  }

  /*
  Graph Bulding Steps (for CSR):
    - Read edgelist once to determine vertex degrees (CountDegrees)
    - Determine vertex offsets by a prefix sum (ParallelPrefixSum)
    - Allocate storage and set points according to offsets (GenIndex)
    - Copy edges into storage
  */
  void MakeCSR(const EdgeList &el, bool transpose, DestID_*** index,
               DestID_** neighs, Partition<NodeID_>* vp, long* pSync, long* pWrk) {
    SGOffset neighbor;
    int64_t receiver, local_v;
    *vp = Partition<NodeID_>(num_nodes_);                                                 // bounds for dividing vertices between PEs
    pvector<NodeID_> degrees = CountDegrees(el, transpose, vp);                  // each pe maintains an array of degrees for the vertices assigned to that pe
    shmem_barrier_all();
    pvector<SGOffset> offsets = ParallelPrefixSum(degrees);                     // offset from start of local neighs array, NOT the global array (symmetric & unsynched)
    //SGOffset* local_max = (SGOffset *) shmem_malloc(sizeof(SGOffset));
    long* local_max = (long *) shmem_malloc(sizeof(long));
    *local_max = (long) *(offsets.begin()+(vp->end-vp->start));
    //SGOffset* max_neigh = (SGOffset *) shmem_malloc(sizeof(SGOffset));
    SGOffset* max_neighbors = (SGOffset *) shmem_calloc(vp->npes, sizeof(SGOffset));    // alternative appproach since long_max_to_all sometimes breaks
    long* max_neigh = (long *) shmem_malloc(sizeof(long));
    //shmem_long_max_to_all(max_neigh, local_max, 1, 0, 0, vp->npes, pWrk, pSync); //all pes must have symmetric neigh arrays, but different lengths. so use max size and acccept the wasted space?               
    for (int i = 0; i < vp->npes; i++)                                          // populate an array of local maxes on every PE
      shmem_putmem(max_neighbors+(vp->pe), local_max, sizeof(SGOffset), i); 
    shmem_barrier_all();
    SGOffset maxn = -1;
    for (int i = 0; i < vp->npes; i++) {
      if (max_neighbors[i] > maxn)
        maxn = max_neighbors[i];
    }
    shmem_barrier_all();
    *neighs = (DestID_ *) shmem_calloc(maxn, sizeof(DestID_));              // maybe copy the neighbors into exact-sized local memory arrays once they no longer need to be symmetric
    //*neighs = (DestID_ *) shmem_calloc(*max_neigh, sizeof(DestID_));              // maybe copy the neighbors into exact-sized local memory arrays once they no longer need to be symmetric
    shmem_barrier_all();
    *index = CSRGraph<NodeID_, DestID_>::GenIndex(offsets, *neighs, vp);
    shmem_barrier_all();
    for (auto it = el.begin(); it < el.end(); it++) {                           // if u || v are part of processing PE's partition, edge must be included on that PE
      Edge e = *it;
      if (symmetrize_ || (!symmetrize_ && !transpose)) {
        receiver = vp->recv(e.u);
        local_v = vp->local_pos(e.u);
        neighbor = shmem_long_atomic_fetch_inc(offsets.begin()+local_v, receiver);                      // move offset pointer for e.u forward one on specified pe
        shmem_putmem((*neighs)+neighbor, &(e.v), sizeof(DestID_), receiver);
      }
      if (symmetrize_ || (!symmetrize_ && transpose)) {
        receiver = vp->recv(e.v);
        local_v = vp->local_pos(static_cast<NodeID_>(e.v));
        // why wont shmem_int64 work?
        neighbor = shmem_long_atomic_fetch_inc(offsets.begin()+local_v, receiver);                     
        DestID_ src = GetSource(e);
        shmem_putmem((*neighs)+neighbor, &src, sizeof(DestID_), receiver);
      }
    }
    shmem_barrier_all();
  }

  CSRGraph<NodeID_, DestID_, invert> MakeGraphFromEL(EdgeList &el, Partition<NodeID_>* p, long* pSync, long* pWrk, int src_opt) {
    DestID_ **index = nullptr, **inv_index = nullptr;
    DestID_ *neighs = nullptr, *inv_neighs = nullptr;
    Timer t;
    t.Start();
    if (num_nodes_ == -1)
      num_nodes_ = FindMaxNodeID(el)+1;
    shmem_barrier_all();
    if (needs_weights_)
      Generator<NodeID_, DestID_, WeightT_>::InsertWeights(el, src_opt);
    shmem_barrier_all();
    MakeCSR(el, false, &index, &neighs, p, pSync, pWrk);
    if (!symmetrize_ && invert) 
      MakeCSR(el, true, &inv_index, &inv_neighs, p, pSync, pWrk);
    shmem_barrier_all();
    t.Stop();
    PrintTime("Build Time", t.Seconds());
    if (symmetrize_)
      return CSRGraph<NodeID_, DestID_, invert>(num_nodes_, index, neighs, pSync, pWrk);
    else
      return CSRGraph<NodeID_, DestID_, invert>(num_nodes_, index, neighs,
                                                inv_index, inv_neighs, pSync, pWrk);
  }

  CSRGraph<NodeID_, DestID_, invert> MakeGraph(long* pWrk, long* pSync) {
    int src_option;                                                             // was the EL from a el file or synthetic? matters for weight generation
    Partition<NodeID_> p; 
    CSRGraph<NodeID_, DestID_, invert> g;
    {  // extra scope to trigger earlier deletion of el (save memory)
      EdgeList el;                                                              // each PE is assigned a portion of the edge list with rr partitioning
      if (cli_.filename() != "") {
        Reader<NodeID_, DestID_, WeightT_, invert> r(cli_.filename());
        if ((r.GetSuffix() == ".sg") || (r.GetSuffix() == ".wsg")) {
          return r.ReadSerializedGraph(pSync, pWrk);
        } else {
          src_option = 0;
          el = r.ReadFile(needs_weights_);
        }
      } else if (cli_.scale() != -1) {
        src_option = 1;
        Generator<NodeID_, DestID_> gen(cli_.scale(), cli_.degree());
        el = gen.GenerateEL(cli_.uniform());
      }
      shmem_barrier_all();
      g = MakeGraphFromEL(el, &p, pSync, pWrk, src_option);
      shmem_barrier_all();
    }
    return SquishGraph(g, &p, pSync, pWrk);
  }

  // distributed k-way merge
  static 
  CSRGraph<NodeID_, DestID_, invert> RelabelByDegree(
      const CSRGraph<NodeID_, DestID_, invert> &g, long* pSync, long* pWrk) {
    if (g.directed()) {
      std::cout << "Cannot relabel directed graph" << std::endl;
      std::exit(-11);
    }
    printf("Rebuilding the graph\n");
    Timer t;
    t.Start();
    // Phase 1: Sort partitioned vectors by degree 
    Partition<NodeID_> vp(g.num_nodes());
    typedef std::pair<int64_t, NodeID_> degree_node_p;
    pvector<degree_node_p> degree_id_pairs(vp.max_width, true);                               // symmetric partitioned array of <node, degree> pairs
    pvector<degree_node_p> temp_pairs(vp.max_width, true);                                    // symmetric partitioned array for storing merged result
    for (NodeID_ n = vp.start; n < vp.end; n++) 
      degree_id_pairs[vp.local_pos(n)] = std::make_pair(g.out_degree(n), n);
    degree_id_pairs.set_widths(vp.max_width, vp.end-vp.start);                                       // Record how many pairs each PE maintains 
    std::sort(degree_id_pairs.begin(), degree_id_pairs.end(), std::greater<degree_node_p>()); // Sort local partition of array
    // Phase 2: K-way merge with tournament trees 
    int* LEADER = (int *) shmem_malloc(sizeof(int));                                       // Each PE is the leader while filling their own temp array
    *LEADER = 0;
    degree_node_p* init_leaves = (degree_node_p *) shmem_calloc(vp.npes, sizeof(degree_node_p));    // Initialize tree with first element from the sorted list on each PE
    if (vp.pe > 0) {
      if (vp.end - vp.start == 0)                                                           // Only PEs with at least one element contribute a non-infinite leaf 
        degree_id_pairs[0] = std::make_pair(std::numeric_limits<int64_t>::max(), 0);
      shmem_putmem(init_leaves+vp.pe, degree_id_pairs.begin(), sizeof(degree_node_p), 0);
    } else {
      init_leaves[0] = degree_id_pairs[0];
    }
    shmem_barrier_all();
    TournamentTree tree(init_leaves, degree_id_pairs);
    shmem_int_wait_until(LEADER, SHMEM_CMP_EQ, vp.pe);                                        // wait until previous PE puts your pe # in LEADER
    for (int i = vp.start; i < vp.end; i++) {                                                  // Leader fills their own temp_pairs before passing leadership to next PE
      temp_pairs[i-vp.start] = tree.pop_root();
      //printf("(%d, %d)\n", temp_pairs[i-vp.start].first, temp_pairs[i-vp.start].second);
    }
    //printf("Pe %d has filled list\n", vp.pe);
    if (vp.pe < vp.npes-1) {
      //printf("PE %d is transfering\n", vp.pe);
      tree.transfer(vp.pe);                                                                    // Transfer contents of tree to next PE
      shmem_int_p(LEADER, vp.pe+1, vp.pe+1);
    }
    //tree.print_tree();
    shmem_barrier_all();
    //shmem_global_exit(0);
    //exit(0);
    degree_id_pairs.~pvector();                                                                         // Free partially sorted lists
    // Phase 3: Relabel vertices by ascending degree 
    pvector<NodeID_> degrees(vp.max_width);
    pvector<NodeID_> new_ids(vp.max_width, true);
    for (NodeID_ n=vp.start; n < vp.end; n++) {
      NodeID_ lp_v = vp.local_pos(n);
      degrees[lp_v] = temp_pairs[lp_v].first;
      shmem_putmem(new_ids.begin()+vp.local_pos(temp_pairs[lp_v].second), &n, sizeof(NodeID_), vp.recv(temp_pairs[lp_v].second));
    }
    shmem_barrier_all();
    //for (int i = vp.start; i < vp.end; i++)
      //printf("PE %d | ID: %d Degree: %d\n", vp.pe, new_ids[vp.local_pos(i)], degrees[vp.local_pos(i)]);
    // Phase 4: Rebuild graph with new IDs 
    pvector<SGOffset> offsets = ParallelPrefixSum(degrees);
    SGOffset* max_neigh = (SGOffset *) shmem_malloc(sizeof(SGOffset));
    shmem_long_max_to_all(max_neigh, offsets.begin()+(vp.end - vp.start), 1, 0, 0, vp.npes, pWrk, pSync);
    DestID_* neighs = (DestID_ *) shmem_calloc(*max_neigh, sizeof(DestID_));
    DestID_** index = CSRGraph<NodeID_, DestID_>::GenIndex(offsets, neighs, &vp);
    shmem_barrier_all();
    for (NodeID_ u = vp.start; u < vp.end; u++) {
      for (NodeID_ v : g.out_neigh(u)) {
        int owner = vp.recv(new_ids[u]);
        int x = new_ids[u] - (owner * vp.partition_width);
        SGOffset off = shmem_long_atomic_fetch_inc(offsets.begin()+x, owner);
        if (owner == vp.pe) {
          neighs[off] = new_ids[v];
        } else {
          shmem_putmem(neighs+off, new_ids.begin()+v, sizeof(DestID_), owner);
        }
      }
    }
    shmem_barrier_all();
    for (int i = 0; i < vp.end-vp.start; i++)
      std::sort(index[i], index[i+1]);
    t.Stop();
    shmem_free(max_neigh);
    PrintTime("Relabel", t.Seconds());
    return CSRGraph<NodeID_, DestID_, invert>(g.num_nodes(), index, neighs, pSync, pWrk);
  }

// UNOPTIMIZED! Possibly requires partitioned parallel sorting and partitioned pvectors for space efficiency.
 /* static
  CSRGraph<NodeID_, DestID_, invert> RelabelByDegree(
      const CSRGraph<NodeID_, DestID_, invert> &g, long* pSync, long* pWrk) {
    shmem_barrier_all();
    if (g.directed()) {
      std::cout << "Cannot relabel directed graph" << std::endl;
      std::exit(-11);
    }
    printf("Rebuilding the graph\n");
    Timer t;
    t.Start();
    Partition<NodeID_> vp(g.num_nodes());
    typedef std::pair<int64_t, NodeID_> degree_node_p;
    pvector<degree_node_p> degree_id_pairs(g.num_nodes(), true);        // complete symmetric list of ids and their degrees, up to date on PE 0

    for (NodeID_ n = vp.start; n < vp.end; n++)
      degree_id_pairs[n] = std::make_pair(g.out_degree(n), n);
    shmem_barrier_all();
    shmem_putmem(degree_id_pairs.begin()+vp.start, degree_id_pairs.begin()+vp.start, (vp.end - vp.start)*sizeof(degree_node_p), 0); // collect on PE 0
    shmem_barrier_all();
    if (vp.pe == 0) {
      std::sort(degree_id_pairs.begin(), degree_id_pairs.end(), std::greater<degree_node_p>());           // sort on PE 0
    //shmem_broadcast32(degree_id_pairs.begin(), degree_id_pairs.begin(), g.num_nodes(), 0, 0, 0, vp.npes, pSync);        // broadcast sorted list
      for (int i = 1; i < vp.npes; i++)
          shmem_putmem(degree_id_pairs.begin(), degree_id_pairs.begin(), g.num_nodes()*sizeof(degree_node_p), i);
    }
    shmem_barrier_all();
    pvector<NodeID_> temp_degrees(g.num_nodes());
    pvector<NodeID_> temp_new_ids(g.num_nodes());
    for (NodeID_ n=0; n < g.num_nodes(); n++) {
      temp_degrees[n] = degree_id_pairs[n].first;
      temp_new_ids[degree_id_pairs[n].second] = n;
    }
    shmem_barrier_all();

    pvector<NodeID_> degrees(vp.max_width);
    pvector<NodeID_> new_ids(vp.max_width);
    for (int i = vp.start; i < vp.end; i++) {
      degrees[i-vp.start] = temp_degrees[i];
      new_ids[i-vp.start] = temp_new_ids[i];
    }

    pvector<SGOffset> offsets = ParallelPrefixSum(degrees);

    SGOffset* max_neigh = (SGOffset *) shmem_malloc(sizeof(SGOffset));
    shmem_long_max_to_all(max_neigh, offsets.begin()+(vp.end - vp.start), 1, 0, 0, vp.npes, pWrk, pSync); 
    DestID_* neighs = (DestID_ *) shmem_calloc(*max_neigh, sizeof(DestID_));
    DestID_** index = CSRGraph<NodeID_, DestID_>::GenIndex(offsets, neighs, &vp);

    shmem_barrier_all();
    for (NodeID_ u = vp.start; u < vp.end; u++) {
      for (NodeID_ v : g.out_neigh(u)) {
        int owner = vp.recv(temp_new_ids[u]);
        int x = temp_new_ids[u] - (owner * vp.partition_width);
        SGOffset off = shmem_long_atomic_fetch_inc(offsets.begin()+x, owner);
        if (owner == vp.pe) {
          neighs[off] = temp_new_ids[v];
        } else {
          shmem_putmem(neighs+off, temp_new_ids.begin()+v, sizeof(DestID_), owner);
        }
      }
    }
    shmem_barrier_all();
    for (int i = 0; i < vp.end-vp.start; i++)
      std::sort(index[i], index[i+1]);

    shmem_barrier_all(); 
    t.Stop();
    shmem_free(max_neigh);
    PrintTime("Relabel", t.Seconds());
    return CSRGraph<NodeID_, DestID_, invert>(g.num_nodes(), index, neighs, pSync, pWrk);
  }*/

};

// complicated optimized partitioned version - https://sites.cs.ucsb.edu/~gilbert/reports/sorting111904.pdf
/*static
CSRGraph<NodeID_, DestID_, invert> RelabelByDegree(
    const CSRGraph<NodeID_, DestID_, invert> &g) {
  if (g.directed()) {
    std::cout << "Cannot relabel directed graph" << std::endl;
    std::exit(-11);
  }
  Timer t;
  t.Start();
  Partition<NodeID_> vp(g.num_nodes());
  typedef std::pair<int64_t, NodeID_> degree_node_p;            // do we need a different partition scheme to satisfy distribution requirement?
  pvector<degree_node_p> degree_id_pairs(vp.max_width);         // symmetric partitioned pvector (unsynched)
  for (NodeID_ n = vp.start; n < vp.end; n++)                   // 1. Local sort
    degree_id_pair(g.out_degree(n), n);                         
    std::sort(degree_id_pairs.begin(), degree_id_pairs.end(),
              std::greater<degree_node_p>());
  shmem_barrier_all();
  // 2. Exact splitting
*/

#endif  // BUILDER_H_
