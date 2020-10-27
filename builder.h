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

  void BAIL() {
    shmem_barrier_all();
    shmem_global_exit(0);
    exit(0);
  }

  DestID_ GetSource(EdgePair<NodeID_, NodeID_> e) {
    return e.u;
  }

  DestID_ GetSource(EdgePair<NodeID_, NodeWeight<NodeID_, WeightT_>> e) {
    return NodeWeight<NodeID_, WeightT_>(e.u, e.v.w);
  }

  NodeID_ FindMaxNodeID(const EdgeList &el) {
    //if (el.symmetric()) {
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
  /*  } else {
      NodeID_ max_seen = 0;
      #pragma omp parallel for reduction(max : max_seen)
      for (auto it = el.begin(); it < el.end(); it++) {
        Edge e = *it;
        max_seen = std::max(max_seen, e.u);
        max_seen = std::max(max_seen, (NodeID_) e.v);
      }
      return max_seen;
    }*/
  }

  // Return pvector representing degrees for vertices assigned to local PE
  // pvector is symmetric and up-to-date but unsynched - do not synch!
  // pvectors on each pe should be concatenated to make a complete list (once unused remainder is trimmed off PEs != npes-1)
  // ASSUMES NODEIDS ARE INTS
  pvector<NodeID_> CountDegrees(const EdgeList &el, bool transpose, Partition<NodeID_>* vp, Partition<>* ep = NULL) {
    int local_v, receiver;
    pvector<NodeID_> degrees(vp->max_width, 0, true);                                     // Symmetric pvector of size max partition width
    printf("Pe %d has max width %d and NodeID_ size %lu\n", vp->pe, vp->max_width, sizeof(NodeID_));
    //#pragma omp parallel for
    Edge e;
    /*shmem_barrier_all();
    shmem_global_exit(0);
    exit(0);*/
    for (auto it = el.begin()/*+(ep->start)*/; it < el.end()/*.begin()+(ep->end)*/; it++) {
      e = *it;
      if (symmetrize_ || (!symmetrize_ && !transpose)) {
        receiver = vp->recv(e.u);
        local_v = vp->local_pos(e.u);
        shmem_int_atomic_inc(degrees.begin()+(local_v), receiver);                                   // increment degree of vertex e.u on pe receiver (could be local PE) 
      }
      if (symmetrize_ || (!symmetrize_ && transpose)) {
        receiver = vp->recv(e.v);
        local_v = vp->local_pos(e.v);
      //  printf("PE %d is about to increment degrees[%d] (node %d) on PE %d\n", vp->pe, local_v, e.v, receiver);
        shmem_int_atomic_inc(degrees.begin()+(local_v), receiver);                                   // increment degree of vertex e.v on pe receiver (could be local PE) 
      }
    }
    printf("I escaped!\n");
    shmem_barrier_all();
    //shmem_global_exit(0);
    //exit(0);
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
  //  printf("PE %d allocating prefix with size %lu\n", shmem_my_pe(), degrees.size());
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
    int indx;
    pvector<NodeID_> diffs(vp.max_width);
    DestID_ *n_start, *n_end;
    //#pragma omp parallel for private(n_start, n_end)
    for (NodeID_ n = vp.start; n < vp.end; n++) {
      indx = n - vp.start;
//      if (vp.pe == 0)
  //      printf("PE %d has indx %d and n %d\n", vp.pe, indx, n);
      if (transpose) {
        n_start = g.in_neigh(n).begin();
        n_end = g.in_neigh(n).end();
      } else {
        n_start = g.out_neigh(n).begin();
        n_end = g.out_neigh(n).end();
      }
      std::sort(n_start, n_end);
      DestID_ *new_end = std::unique(n_start, n_end);
      new_end = std::remove(n_start, new_end, n);
      diffs[indx] = new_end - n_start;
      //if (vp.pe == 0)
        //printf("PE %d has new_end %d and n_start %d\n", vp.pe, *new_end, *n_start);
    }
    //printf("PE %d says diffs.size = %lu\n", vp.pe, diffs.size());
    pvector<SGOffset> sq_offsets = ParallelPrefixSum(diffs);
    SGOffset* max_neigh = (SGOffset *) shmem_malloc(sizeof(SGOffset));
    shmem_long_max_to_all(max_neigh, sq_offsets.begin()+(vp.end - vp.start), 1, 0, 0, vp.npes, pWrk, pSync); 
    //*sq_neighs = new DestID_[sq_offsets[vp.end-vp.start]];
    //printf("pe %d - squish neigh len = %lu\n", vp.pe, sq_offsets[vp.end-vp.start]);
    *sq_neighs = (DestID_ *) shmem_calloc(*max_neigh, sizeof(DestID_));
    *sq_index = CSRGraph<NodeID_, DestID_>::GenIndex(sq_offsets, *sq_neighs, &vp);
    //printf("PE %d | fking sq_indx: %p\n", vp.pe, (void*) *sq_index);
    #pragma omp parallel for private(n_start)
    for (NodeID_ n=vp.start; n < vp.end; n++) {
      indx = n - vp.start;
      if (transpose)
        n_start = g.in_neigh(n).begin();
      else
        n_start = g.out_neigh(n).begin();
      std::copy(n_start, n_start+diffs[indx], (*sq_index)[indx]);
    }
  }

  CSRGraph<NodeID_, DestID_, invert> SquishGraph(
      const CSRGraph<NodeID_, DestID_, invert> &g, Partition<NodeID_>* vp, long* pSync, long* pWrk) {
    DestID_ **out_index, *out_neighs, **in_index, *in_neighs;
    SquishCSR(g, false, &out_index, &out_neighs, *vp, pSync, pWrk);
    //shmem_barrier_all();
    //shmem_set_lock(PRINT_LOCK);
    //for (int i = vp->start; i <= vp->end; i++)
      //printf("PE: %d Node %d has out_index %p\n", vp->pe, i, (void *) out_index[i-vp->start]); /* *(out_index[i-vp->start]));*/
    //shmem_clear_lock(PRINT_LOCK);
    //printf("PE %d has out_index start: %p\n", vp->pe, (void*) &out_index);
    if (g.directed()) {
      if (invert)
        SquishCSR(g, true, &in_index, &in_neighs, *vp, pSync, pWrk);
      //printf("PE %d | out index: %p => %d | in index: %p => %d\n", vp->pe, (void *) out_index, **out_index, (void *) in_index, **in_index);
      return CSRGraph<NodeID_, DestID_, invert>(g.num_nodes(), out_index,
                                                out_neighs, in_index,
                                                in_neighs, pSync, pWrk);
    } else {
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
    printf("Entered MakeCSR\n");
    SGOffset neighbor;
    int64_t receiver, local_v;
    *vp = Partition<NodeID_>(num_nodes_);                                                 // bounds for dividing vertices between PEs
    printf("PE %d has partition start = %d and end = %d\n", vp->pe, vp->start, vp->end);
    //long* num_edges = (long*) shmem_calloc(1, sizeof(long));
    //shmem_long_sum_to_all(num_edges, el.local_width(), 1, 0, 0, shmem_n_pes(), pWrk, pSync);  
    //Partition<> ep(el.size());                                                     // bounds for dividing work (edges to process) between PEs
    pvector<NodeID_> degrees = CountDegrees(el, transpose, vp);                  // each pe maintains an array of degrees for the vertices assigned to that pe
    printf("Degrees size = %lu\n", degrees.size());
    //shmem_global_exit(0);
    //exit(0);
    shmem_barrier_all();
    pvector<SGOffset> offsets = ParallelPrefixSum(degrees);                     // offset from start of local neighs array, NOT the global array (symmetric & unsynched)
    /*if (transpose) {
      int i = 0;
      for (d : degrees) {
        printf("PE %d | Node %d has degree %lu and offset %lu\n", vp->pe, (vp->start)+i, degrees[i], offsets[i]);
        i++;
      }
      printf("PE %d | Final offset: %lu\n", vp->pe, offsets[i++]);
      shmem_barrier_all();
      //BAIL();
    }*/
    //printf("Offsets.size(): %lu, final offset: %lu\n", offsets.size(), offsets[vp->end - vp->start]);
    SGOffset* max_neigh = (SGOffset *) shmem_malloc(sizeof(SGOffset));
    shmem_long_max_to_all(max_neigh, offsets.begin()+(vp->end - vp->start), 1, 0, 0, vp->npes, pWrk, pSync); //all pes must have symmetric neigh arrays, but different lengths. so use max size and acccept the wasted space?               

    printf("oh boy here i go killin again\n");
    //shmem_global_exit(0);
    //exit(0);
    *neighs = (DestID_ *) shmem_calloc(*max_neigh, sizeof(DestID_));              // maybe copy the neighbors into exact-sized local memory arrays once they no longer need to be symmetric
    *index = CSRGraph<NodeID_, DestID_>::GenIndex(offsets, *neighs, vp);
    shmem_barrier_all();
    for (auto it = el.begin()/*+ep.start*/; it < el.end()/*el.begin()+ep.end*/; it++) {                           // if u || v are part of processing PE's partition, edge must be included on that PE
      Edge e = *it;
      if (symmetrize_ || (!symmetrize_ && !transpose)) {
        receiver = vp->recv(e.u);
        local_v = vp->local_pos(e.u);
        neighbor = shmem_long_atomic_fetch_inc(offsets.begin()+local_v, receiver);                      // move offset pointer for e.u forward one on specified pe
        //shmem_int_put((*neighs)+neighbor, &(e.v), 1, receiver);                        // should be safe since offset pointer has already been atomically incremented?
        shmem_putmem((*neighs)+neighbor, &(e.v), sizeof(DestID_), receiver);
          //(*neighs)[fetch_and_add(offsets[e.u], 1)] = e.v;                        // array of neighbors[old offset for e.u] = e.v, increment old offset for e.u
      }
      //shmem_barrier_all();
      if (symmetrize_ || (!symmetrize_ && transpose)) {
        receiver = vp->recv(e.v);
        local_v = vp->local_pos(static_cast<NodeID_>(e.v));
        // why wont shmem_int64 work?
        neighbor = shmem_long_atomic_fetch_inc(offsets.begin()+local_v, receiver);                     
        NodeID_ src = GetSource(e);
        //shmem_int_put((*neighs)+neighbor, &src, 1, receiver);                        
        shmem_putmem((*neighs)+neighbor, &src, sizeof(DestID_), receiver);
      }
      //shmem_barrier_all();              // are these necessary?
    }


// when/how best to free the symmetric neighbor array?

    //printf("PE %d: (%p => %d)\n", vp->pe, (void*) **index, ***index);
    //printf("PE %d: (%p => %d)\n", vp->pe, (void*) *((*index)+1), **((*index)+1));
    /*shmem_barrier_all();
    if (transpose) {
      int i = 0;
      for (auto it = vp->start; it < vp->end; it++) {
        printf("PE %d | Node %d has neighbors ", vp->pe, vp->global_pos(i));
        for (o : offsets)
          printf("%d ", *((*neighs)+o));
        i++;
      }
    }*/
    /*printf("PE %d | ", vp->pe);
    for (int i = 0; i < 7; i++)
      printf("%d ", (*neighs)[i]);
    printf("\n");
    printf("Coheck 2\n");*/
  }

  CSRGraph<NodeID_, DestID_, invert> MakeGraphFromEL(EdgeList &el, Partition<NodeID_>* p, long* pSync, long* pWrk, int src_opt) {
    DestID_ **index = nullptr, **inv_index = nullptr;
    DestID_ *neighs = nullptr, *inv_neighs = nullptr;
    Timer t;
    t.Start();
    if (num_nodes_ == -1)
      num_nodes_ = FindMaxNodeID(el)+1;
    shmem_barrier_all();
    printf("Max node is %lu\n", num_nodes_);
    if (needs_weights_) {
      printf("ing weights\n");
      Generator<NodeID_, DestID_, WeightT_>::InsertWeights(el, src_opt);
    }
    printf("PE %d | EL size inside MakeGraphFromEL: %lu\n", shmem_my_pe(), el.size());
    shmem_barrier_all();
    // transplanted logic
    int local_v, receiver;
    *p = Partition<NodeID_>(num_nodes_);                                                 // bounds for dividing vertices between PEs
    pvector<NodeID_> degrees(p->max_width, 0, true);                                     // Symmetric pvector of size max partition width
    printf("Pe %d has max width %d and NodeID_ size %lu\n", p->pe, p->max_width, sizeof(NodeID_));
    Edge e;
    for (auto it = el.begin()/*+(ep->start)*/; it < el.end()/*.begin()+(ep->end)*/; it++) {
      e = *it;
      if (symmetrize_ || (!symmetrize_ && !false)) {
        receiver = p->recv(e.u);
        local_v = p->local_pos(e.u);
        shmem_int_atomic_inc(degrees.begin()+(local_v), receiver);                                   // increment degree of vertex e.u on pe receiver (could be local PE) 
      }
      if (symmetrize_ || (!symmetrize_ && false)) {
        receiver = p->recv(e.v);
        local_v = p->local_pos(e.v);
      //  printf("PE %d is about to increment degrees[%d] (node %d) on PE %d\n", vp->pe, local_v, e.v, receiver);
        shmem_int_atomic_inc(degrees.begin()+(local_v), receiver);                                   // increment degree of vertex e.v on pe receiver (could be local PE) 
      }
    }
    printf("I escaped!\n");
    BAIL();
    // end transplant
   // shmem_global_exit(0);
    //exit(0);
    MakeCSR(el, false, &index, &neighs, p, pSync, pWrk);
    printf("Made CSR\n");
  //  shmem_global_exit(0);
  //  exit(0);
    if (!symmetrize_ && invert) {
      MakeCSR(el, true, &inv_index, &inv_neighs, p, pSync, pWrk);
    }
    //printf("check\n");
    //std::cout << std::flush;
    shmem_barrier_all();
    //BAIL();
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
          printf("more are coming\n");
        //  ~el;
        //  shmem_global_exit(0);
        //  exit(0);
        }
      } else if (cli_.scale() != -1) {
        src_option = 1;
        Generator<NodeID_, DestID_> gen(cli_.scale(), cli_.degree());
        el = gen.GenerateEL(cli_.uniform());
        printf("check\n");
        //BAIL();
        shmem_barrier_all();
      }
    //  shmem_barrier_all();
     // printf("Here\n");
    //} // extra scope, remove later
    //shmem_barrier_all();
    //shmem_global_exit(0);
    //exit(0);
      shmem_barrier_all();
      printf("EL is created with size %lu, EL[0]: (%d, %d)\n", el.size(), (*(el.begin())).u, (*(el.begin())).v);
      //for (EdgePair<int, int> x : el)
      //  printf("PE %d | (%d, %d)\n", shmem_my_pe(), x.u, x.v);
      //shmem_barrier_all();
      g = MakeGraphFromEL(el, &p, pSync, pWrk, src_option);
    }
    shmem_barrier_all();
    shmem_global_exit(0);
    exit(0);
    return SquishGraph(g, &p, pSync, pWrk);
  }

// UNOPTIMIZED! Possibly requires partitioned parallel sorting and partitioned pvectors for space efficiency.
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
    Partition<NodeID_> vp(g.num_nodes());
    typedef std::pair<int64_t, NodeID_> degree_node_p;
    pvector<degree_node_p> degree_id_pairs(g.num_nodes(), true);        // up to date on PE 0
    for (NodeID_ n = vp.start; n < vp.end; n++)
      degree_id_pairs[n] = std::make_pair(g.out_degree(n), n);
    shmem_putmem(degree_id_pairs.begin()+vp.start, degree_id_pairs.begin()+vp.start, (vp.end - vp.start)*sizeof(degree_node_p), 0); // collect on PE 0
    shmem_barrier_all();
    if (vp.pe == 0) {
      std::sort(degree_id_pairs.begin(), degree_id_pairs.end(), std::greater<degree_node_p>());           // sort on PE 0
    //shmem_broadcast32(degree_id_pairs.begin(), degree_id_pairs.begin(), g.num_nodes(), 0, 0, 0, vp.npes, pSync);        // broadcast sorted list
      for (int i = 0; i < vp.npes; i++)
        if (vp.pe != i)
          shmem_putmem(degree_id_pairs.begin(), degree_id_pairs.begin(), g.num_nodes()*sizeof(degree_node_p), i);
    }
    shmem_barrier_all();
    //if (vp.pe == 1)
      //for (auto it : degree_id_pairs)
        //printf("degree(%d) = %lu\n", it.second, it.first);
    printf("Check x\n");
    pvector<NodeID_> temp_degrees(g.num_nodes());
    pvector<NodeID_> temp_new_ids(g.num_nodes());
    for (NodeID_ n=0; n < g.num_nodes(); n++) {
      temp_degrees[n] = degree_id_pairs[n].first;
      temp_new_ids[degree_id_pairs[n].second] = n;
    }
    shmem_barrier_all();
   // for (int i = 0; i < g.num_nodes(); i++)
   //   printf("PE %d | degree = %d, NodeID = %d\n", vp.pe, temp
    pvector<NodeID_> degrees(vp.max_width);
    pvector<NodeID_> new_ids(vp.max_width);
    for (int i = vp.start; i < vp.end; i++) {
      degrees[i-vp.start] = temp_degrees[i];
      new_ids[i-vp.start] = temp_new_ids[i];
    }
    //for (int i = 0; i < vp.max_width; i++)
      //printf("PE %d | degrees[%d] = %d, new_ids[%d] = %d\n", vp.pe, i, degrees[i], i, new_ids[i]);
    pvector<SGOffset> offsets = ParallelPrefixSum(degrees);
    //for (o : offsets)
      //printf("Pe: %d | o = %d\n", vp.pe, o);
    SGOffset* max_neigh = (SGOffset *) shmem_malloc(sizeof(SGOffset));
    shmem_long_max_to_all(max_neigh, offsets.begin()+(vp.end - vp.start), 1, 0, 0, vp.npes, pWrk, pSync); 
    printf("PE %d | Max neigh: %lu\n", vp.pe, *max_neigh);
    DestID_* neighs = (DestID_ *) shmem_calloc(*max_neigh, sizeof(DestID_));
    DestID_** index = CSRGraph<NodeID_, DestID_>::GenIndex(offsets, neighs, &vp);
    shmem_barrier_all();
    for (NodeID_ u = vp.start; u < vp.end; u++) {
      for (NodeID_ v : g.out_neigh(u)) {
        //printf("PE %d | u = %d, v = %d, new_ids[u] = %d\n", vp.pe, u, v, temp_new_ids[u]);
        //if (temp_new_ids[u]
        //if (temp_new_ids[u]   
        //neighs[offsets[temp_new_ids[u]-vp.start]+k] = temp_new_ids[v];
        //if (temp_new_ids[u] >= vp.start && temp_new_ids < vp.end) 
        int owner = vp.recv(temp_new_ids[u]);
        int x = temp_new_ids[u] - (owner * vp.partition_width);
        SGOffset off = shmem_long_atomic_fetch_inc(offsets.begin()+x, owner);
        //printf("Pe %d | U = %d | V = %d | owner = %d | offset = %lu\n", vp.pe, u, v, owner, off);
        if (owner == vp.pe) {
      //    printf("U = %d, V = %d | PE %d is locally putting %d at neighs[%d]\n", u, v, vp.pe, temp_new_ids[v], off);
          neighs[off] = temp_new_ids[v];
        } else {
      //    printf("U = %d, V = %d | PE %d is remotely putting %d at neighs[%d] on PE %d\n", u, v, vp.pe, temp_new_ids[v], off, owner);
          shmem_putmem(neighs+off, temp_new_ids.begin()+v, sizeof(DestID_), owner);
        }
        //printf("PE %d | u = %d, v = %d, new_id(u) = %d, new_id(v) = %d, offsets[%d] = %d\n", vp.pe, u, v, temp_new_ids[u], temp_new_ids[v], temp_new_ids[u], off);
        //  k++;
        //printf("Current offset of %d = %lu\n", temp_new_ids[u], offsets[(temp_new_ids[u]-vp.start)]);
        //SGOffset off = offsets
      //std::sort(index[new_ids[u-vp.start]], index[new_ids[u-vp.start]+1]);
      }
    }
    //printf("hey\n");
    //shmem_clear_lock(LOCK);
    shmem_barrier_all();
    //for (int n = 0; n < *max_neigh; n++)
    //  printf("PE %d | n = %d\n", 
    for (int i = 0; i < vp.end-vp.start/*)+1vp.max_width*/; i++)
      std::sort(index[i], index[i+1]);
    shmem_barrier_all(); 
    //printf("what\n");
    //for (int n = 0; n < *max_neigh; n++)
      //printf("Pe %d | %d\n", vp.pe, neighs[n]);
    t.Stop();
    shmem_free(max_neigh);
    PrintTime("Relabel", t.Seconds());
    //return CSRGraph<NodeID_, DestID_, invert>();
    return CSRGraph<NodeID_, DestID_, invert>(g.num_nodes(), index, neighs, pSync, pWrk);
  }
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
