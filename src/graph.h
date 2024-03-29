// Copyright (c) 2015, The Regents of the University of California (Regents)
// See LICENSE.txt for license details

#ifndef GRAPH_H_
#define GRAPH_H_

#include <algorithm>
#include <cinttypes>
#include <cstddef>
#include <iostream>
#include <type_traits>

#include "pvector.h"
#include "util.h"


/*
GAP Benchmark Suite
Class:  CSRGraph
Author: Scott Beamer

Simple container for graph in CSR format
 - Intended to be constructed by a Builder
 - To make weighted, set DestID_ template type to NodeWeight
 - MakeInverse parameter controls whether graph stores its inverse
*/

/*
DMM-GAPBS
Author: Zach Hansen
Adaptation Notes:
 - Neighborhood methods are overwritten to enable PEs to iterate over neighborhoods
   on foreign PEs using range based iteration (performs shmem gets behind the scenes)
 - shmem ptrs should be used when two PEs share a memory space, but there is no check for that yet
 - Enables printing of the distributed topology 
*/

// Used to hold node & weight, with another node it makes a weighted edge
template <typename NodeID_, typename WeightT_>
struct NodeWeight {
  NodeID_ v;
  WeightT_ w;
  NodeWeight() {}
  NodeWeight(NodeID_ v) : v(v), w(1) {}
  NodeWeight(NodeID_ v, WeightT_ w) : v(v), w(w) {}

  bool operator< (const NodeWeight& rhs) const {
    return v == rhs.v ? w < rhs.w : v < rhs.v;
  }

  // doesn't check WeightT_s, needed to remove duplicate edges
  bool operator== (const NodeWeight& rhs) const {
    return v == rhs.v;
  }

  // doesn't check WeightT_s, needed to remove self edges
  bool operator== (const NodeID_& rhs) const {
    return v == rhs;
  }

  operator NodeID_() {
    return v;
  }
};

template <typename NodeID_, typename WeightT_>
std::ostream& operator<<(std::ostream& os,
                         const NodeWeight<NodeID_, WeightT_>& nw) {
  os << nw.v << " " << nw.w;
  return os;
}

template <typename NodeID_, typename WeightT_>
std::istream& operator>>(std::istream& is, NodeWeight<NodeID_, WeightT_>& nw) {
  is >> nw.v >> nw.w;
  return is;
}



// Syntatic sugar for an edge
template <typename SrcT, typename DstT = SrcT>
struct EdgePair {
  SrcT u;
  DstT v;

  EdgePair() {}

  EdgePair(SrcT u, DstT v) : u(u), v(v) {}
};

// SG = serialized graph, these types are for writing graph to file
typedef int32_t SGID;
typedef EdgePair<SGID> SGEdge;
typedef int64_t SGOffset;

// If 2 PEs share a memory space, it is more efficient to use the shmem_ptrs
// But PEs on separate machines need to use the get/put mems

template <class NodeID_, class DestID_ = NodeID_, bool MakeInverse = true>
class CSRGraph {
  // Used for *non-negative* offsets within a neighborhood
  typedef std::make_unsigned<std::ptrdiff_t>::type OffsetT;

  // Used to access neighbors of vertex, basically sugar for iterators
  class Neighborhood {
    int owner;
    bool same_space_;
    bool bulk_;
    NodeID_ n_;                 
    NodeID_ local;
    DestID_** g_index_;
    DestID_** ptr_start;
    DestID_** ptr_end;
    OffsetT start_offset_;
    DestID_* foreign_start;
    DestID_* foreign_end;
    DestID_* beginning;
    DestID_* ending;
    DestID_* buff_neighbor_;
   public:
    DestID_* current;           // Current position of iterator within neighborhood
    DestID_ curr_val;           // Current value pointed at by iterator

    Neighborhood(NodeID_ n, DestID_** g_index, Partition<NodeID_> vp, OffsetT start_offset, DestID_* buffer, bool bulk_transfer, bool same_space = false) :
        n_(n), g_index_(g_index), start_offset_(0), buff_neighbor_(buffer), bulk_(bulk_transfer), same_space_(same_space) {
      OffsetT max_offset;
      local = vp.local_pos(n_);
      owner = vp.recv(n_);
      //printf("PE %d | buffer: %p\n", shmem_my_pe(), (void*) buff_neighbor_);
      //printf("PE %d | g_index: %p, globeql %d local %d\n", shmem_my_pe(), (void*) g_index_, n, local);
      if (vp.pe == owner) {
        max_offset = g_index_[local+1] - g_index_[local];
        start_offset_ = std::min(start_offset, max_offset);
        beginning = g_index_[local] + start_offset_;
        ending = g_index_[local+1];
        current = beginning;
      } else {
        if (same_space_) {    // true if calling PE shares a memory space with owner pe
          ptr_start = (DestID_**) shmem_ptr(g_index_+local, owner);
          ptr_end = (DestID_ **) shmem_ptr((g_index_+(local+1)), owner);
          max_offset = ptr_end - ptr_start;
          start_offset_ = std::min(start_offset, max_offset);
          beginning = *ptr_start + start_offset_;
          ending = *ptr_end;
          current = beginning;
        } else {
          shmem_getmem(&foreign_start, g_index_+local, sizeof(DestID_*), owner);
          shmem_getmem(&foreign_end, g_index_+(local+1), sizeof(DestID_*), owner);
          max_offset = foreign_end - foreign_start; 
          start_offset_ = std::min(start_offset, max_offset);
          beginning = foreign_start + start_offset_;
          ending = foreign_end;
	  shmem_getmem(buff_neighbor_, beginning, (ending-beginning)*sizeof(DestID_), owner);
	  //printf("HHH: %lu\n", ending - beginning);
          current = buff_neighbor_;
	  ending = buff_neighbor_ + (ending-beginning); 
	  //for (int i = 0; i < ending-beginning; i++) {
//		printf("Neigh: %d\n", current[i]);	  
//	  }
          //printf("PE %d | Node: %d | (beginning = %p) => %d | (ending = %p) => %d\n", vp.pe, n, (void*) (beginning), *beginning, (void*) ending, *ending);
        }
      }
      //current = beginning;
    }

    // begin and end are used for range based iteration, must return refs to neighborhoods instead of DestIDs
    // otherwise the overloaded operators won't work
    Neighborhood& begin() { return *this; }

    Neighborhood& end() { return *this; }

    typedef DestID_* iterator;

    // start and finish are what begin and end used to be: local memory can be directly dereferenced with these,
    // but for accessing PEs on separate computers they can only be used to get the address for use in a different shmem call
    iterator start(bool help = false) { 
      if (shmem_my_pe() == owner) {
        //printf("PE %d is starting for node %d\n", shmem_my_pe(), n_);
        return(beginning); 
      } else {
        if (help)
          printf("PE %d is calling n_start on a foreign PE for node %d\n", shmem_my_pe(), n_);
        if (same_space_) {            // check to see if the PEs share a memory space                                                                                   
          DestID_* neigh_start = (DestID_ *) shmem_ptr(beginning, owner);
          return(neigh_start);
        } else {
          //printf("PE %d is requesting start address %p\n", shmem_my_pe(), (void*) (beginning));
          return(beginning);
        }
      }
    }

    iterator finish(bool help = false) { 
      if (shmem_my_pe() == owner)
        return(ending); 
      else {
        if (same_space_)
          return ((DestID_*) shmem_ptr(ending, owner));
        else
          return(foreign_end);
      }
    }

    bool operator!=(Neighborhood const& it) const { return it.ending != current; }

    DestID_& operator*() {
        return(*current);
    }

    const DestID_& operator*() const {
        return(*current);
    }

    Neighborhood& operator+(size_t n) {
      current = current + n;                    // Is this correct pointer arithmetic?
      return *this;
    }

    // prefix increment operator
    Neighborhood& operator++() {
      ++current;
      return *this;
    }

    // postfix increment operator
    Neighborhood operator++(int) {
      ++current;
      return *this;
    }

    DestID_& operator[](size_t n) {
      if (shmem_my_pe() == owner) {
        return(*(beginning+n));
      } else {
	      if (bulk_) {
		      	return(*(buff_neighbor_+n));
		} else {
        		shmem_getmem(&curr_val, beginning+n, sizeof(DestID_), owner); 
        		return curr_val;
		}
      }
    }

    const DestID_& operator[](size_t n) const {
      if (shmem_my_pe() == owner) {
        return(*(beginning+n));
      } else {
	      if (bulk_) {
		      	return(*(buff_neighbor_+n));
		} else {
        		shmem_getmem(&curr_val, beginning+n, sizeof(DestID_), owner);         
        		return curr_val;
		}
      }
    }
  };

  void ReleaseResources() {
    if (out_index_ != nullptr)
      shmem_free(out_index_);
    if (out_neighbors_ != nullptr)
      shmem_free(out_neighbors_);
    if (directed_) {
      if (in_index_ != nullptr)
        shmem_free(in_index_);
      if (in_neighbors_ != nullptr)
        shmem_free(in_neighbors_);
      if (buff_neighbor_ != nullptr)
        delete [] buff_neighbor_;
    }
  }


 public:
  CSRGraph() : directed_(false), num_nodes_(-1), num_edges_(-1),
    out_index_(nullptr), out_neighbors_(nullptr),
    in_index_(nullptr), in_neighbors_(nullptr) {}

  CSRGraph(int64_t num_nodes, DestID_** index, DestID_* neighs, long* pSync, long* pWrk) :
    directed_(false), num_nodes_(num_nodes),
    out_index_(index), out_neighbors_(neighs),
    in_index_(index), in_neighbors_(neighs) {
    int64_t* edge_counts = (int64_t*) shmem_malloc(sizeof(int64_t));
    Partition<NodeID_> p(num_nodes_);
    *edge_counts = out_index_[p.end - p.start] - out_index_[0];                             // how long is the local neighbor array?
    shmem_long_sum_to_all(edge_counts, edge_counts, 1, 0, 0, p.npes, pWrk, pSync);      // Reduction : +
    num_edges_ = *edge_counts / 2;
    buff_neighbor_ = new DestID_[num_edges_/num_nodes_];
  }

  CSRGraph(int64_t num_nodes, DestID_** out_index, DestID_* out_neighs,
        DestID_** in_index, DestID_* in_neighs, long* pSync, long* pWrk, bool symmetric = true) :
    directed_(true), num_nodes_(num_nodes),
    out_index_(out_index), out_neighbors_(out_neighs),
    in_index_(in_index), in_neighbors_(in_neighs)/*, p_{num_nodes}*/ {
    int64_t* edge_counts = (int64_t*) shmem_malloc(sizeof(int64_t));
    Partition<NodeID_> p(num_nodes_);
    *edge_counts = out_index_[p.end - p.start] - out_index_[0];                             // how long is the local neighbor array?
    shmem_long_sum_to_all(edge_counts, edge_counts, 1, 0, 0, p.npes, pWrk, pSync);      // Reduction : +
    num_edges_ = *edge_counts;
    buff_neighbor_ = new DestID_[num_edges_/num_nodes_];				// is this necessary? can it be inherited from the graph pre-inversion?
  }

  CSRGraph(CSRGraph&& other) : directed_(other.directed_),
    num_nodes_(other.num_nodes_), num_edges_(other.num_edges_),
    out_index_(other.out_index_), out_neighbors_(other.out_neighbors_),
    in_index_(other.in_index_), in_neighbors_(other.in_neighbors_), buff_neighbor_(other.buff_neighbor_) {
      other.num_edges_ = -1;
      other.num_nodes_ = -1;
      other.out_index_ = nullptr;
      other.out_neighbors_ = nullptr;
      other.in_index_ = nullptr;
      other.in_neighbors_ = nullptr;
      other.buff_neighbor_ = nullptr;
  }

  ~CSRGraph() {
    ReleaseResources();
  }

  CSRGraph& operator=(CSRGraph&& other) {
    if (this != &other) {
      ReleaseResources();
      directed_ = other.directed_;
      num_edges_ = other.num_edges_;
      num_nodes_ = other.num_nodes_;
      out_index_ = other.out_index_;
      out_neighbors_ = other.out_neighbors_;
      in_index_ = other.in_index_;
      in_neighbors_ = other.in_neighbors_;
      buff_neighbor_ = other.buff_neighbor_;
      other.num_edges_ = -1;
      other.num_nodes_ = -1;
      other.out_index_ = nullptr;
      other.out_neighbors_ = nullptr;
      other.in_index_ = nullptr;
      other.in_neighbors_ = nullptr;
      other.buff_neighbor_ = nullptr;
    }
    return *this;
  }

  bool directed() const {
    return directed_;
  }

  int64_t num_nodes() const {
    return num_nodes_;
  }

  int64_t num_edges() const {
    return num_edges_;
  }

  int64_t num_edges_directed() const {
    return directed_ ? num_edges_ : 2*num_edges_;
  }

  int64_t out_degree(NodeID_ v, bool same_space = false) const {                // Do the owner of v and the calling PE share a memory space?
    int64_t degree;                                                                                             
    Partition<NodeID_> vp(num_nodes_);
    NodeID_ local = vp.local_pos(v);    
    if (v >= vp.start && v < vp.end) {
      degree = out_index_[local+1] - out_index_[local];
    } else {
      if (same_space) {
        NodeID_** p_one = (NodeID_**) shmem_ptr(out_index_+local, vp.recv(v));
        NodeID_** p_two = (NodeID_ **) shmem_ptr((out_index_+(local+1)), vp.recv(v));
        degree = *p_two - *p_one;
      } else {
        DestID_* one;
        DestID_* two;
        shmem_getmem(&one, out_index_+local, sizeof(DestID_*), vp.recv(v));  
        shmem_getmem(&two, out_index_+(local+1), sizeof(DestID_*), vp.recv(v)); 
        degree = two - one;
      }
    }
    return degree;
  }

  int64_t in_degree(NodeID_ v, bool same_space = false) const {
    static_assert(MakeInverse, "Graph inversion disabled but reading inverse");
    int64_t degree;
    Partition<NodeID_> vp(num_nodes_);
    NodeID_ local = vp.local_pos(v);
    if (v >= vp.start && v < vp.end) {
      degree = in_index_[local+1] - in_index_[local];
    } else {
      if (same_space) {
        NodeID_** p_one = (NodeID_**) shmem_ptr(in_index_+local, vp.recv(v));
        NodeID_** p_two = (NodeID_ **) shmem_ptr((in_index_+(local+1)), vp.recv(v));
        degree = *p_two - *p_one;
      } else {
        DestID_* one;
        DestID_* two;
        shmem_getmem(&one, in_index_+local, sizeof(DestID_*), vp.recv(v));  
        shmem_getmem(&two, in_index_+(local+1), sizeof(DestID_*), vp.recv(v));  
        degree = two - one;
      }
    }
    return degree;
  }

  Neighborhood out_neigh(NodeID_ n, OffsetT start_offset = 0) const {
    Partition<NodeID_> vp(num_nodes_);
    return Neighborhood(n, out_index_, vp, start_offset, buff_neighbor_, true);
  }

  Neighborhood in_neigh(NodeID_ n, OffsetT start_offset = 0) const {
    Partition<NodeID_> vp(num_nodes_);
    static_assert(MakeInverse, "Graph inversion disabled but reading inverse");
    return Neighborhood(n, in_index_, vp, start_offset, buff_neighbor_, true);
  }

  void PrintStats() const {
    std::cout << "Graph has " << num_nodes_ << " nodes and "
              << num_edges_ << " ";
    if (!directed_)
      std::cout << "un";
    std::cout << "directed edges for degree: ";
    std::cout << num_edges_/num_nodes_ << std::endl;
  }

  void PrintTopology(long* PRINT_LOCK) {
    shmem_barrier_all();
    Partition<NodeID_> vp(num_nodes_);
    shmem_set_lock(PRINT_LOCK);
    std::cout << "########################  Graph Topology (Outgoing): PE " << vp.pe <<  " #######################" << std::endl;
    for (NodeID_ i = vp.start; i < vp.end; i++) {
      std::cout << i << ": ";
      //printf("Node %d begin = %d, end = %d\n", i, *(out_neigh(i).begin()), *(out_neigh(i).end()));
      for (DestID_ j : out_neigh(i))
        std::cout << j << " ";
      std::cout << std::endl;
    }
    shmem_clear_lock(PRINT_LOCK);
    shmem_barrier_all();
  }


  void PrintTopology(bool outgoing = true) const {
    Partition<NodeID_> vp(num_nodes_);
    int* PRINTER = (int *) shmem_calloc(1, sizeof(int));                // init 0
    shmem_int_wait_until(PRINTER, SHMEM_CMP_EQ, vp.pe);           // wait until previous PE puts your pe # in PRINTER
    if (outgoing) {
      std::cout << "########################  Graph Topology (Outgoing): PE " << vp.pe <<  " #######################" << std::endl;
      int j = 0;
      for (NodeID_ i = vp.start; i < vp.end; i++) {
        int k = 0;
        //printf("PE %d | i %d | start: %p => %d | end: %p => %d\n", vp.pe, i, (void*) out_neigh(i).start(), out_neigh(i).start()[0], (void*) out_neigh(i).finish(), out_neigh(i).finish()[0]);
        std::cout << i << ": ";
        for (DestID_ j : out_neigh(i)) {
          std::cout << j << " ";
          if (k > 20)
            break;
          k++;
        }
        //std::cout << std::endl;
        printf("\n");
        if (j > 20)
          break;
        j++;
      }
      if (!(vp.pe == vp.npes-1))
        shmem_int_p(PRINTER, vp.pe+1, vp.pe+1);             // who's next?
      shmem_barrier_all();
    } else {
      std::cout << "########################  Graph Topology (Incoming): PE " << vp.pe <<  " #######################" << std::endl;
      int j = 0;
      for (NodeID_ i = vp.start; i < vp.end; i++) {
        int k = 0;
        std::cout << i << ": ";
        for (DestID_ j : in_neigh(i)) {
          std::cout << j << " ";
          if (k > 20)
            break;
          k++;
        }
        //std::cout << std::endl;
        printf("\n");
        if (j > 20)
          break;
        j++;
      }
      if (!(vp.pe == vp.npes-1))
        shmem_int_p(PRINTER, vp.pe+1, vp.pe+1);             // who's next?
      shmem_barrier_all();
    }
    shmem_free(PRINTER);
  }


  // offsets for given pe start from the pes first elem
  // some unused space - max_width = partition_width + remainder
  static DestID_** GenIndex(const pvector<SGOffset> &offsets, DestID_* neighs, Partition<NodeID_>* p) {
    DestID_** index = (DestID_**) shmem_calloc((p->max_width)+1, sizeof(DestID_*)); 
    //#pragma omp parallel for
    for (NodeID_ n = p->start; n <= p->end; n++) {
      index[n-p->start] = neighs + offsets[n-p->start];
    }
    return index;
  }

  static DestID_** OLDGenIndex(const pvector<SGOffset> &offsets, DestID_* neighs) {
    NodeID_ length = offsets.size();
    DestID_** index = new DestID_*[length];
    #pragma omp parallel for
    for (NodeID_ n=0; n < length; n++)
      index[n] = neighs + offsets[n];
    return index;
  }

  pvector<SGOffset> VertexOffsets(Partition<NodeID_>* vp, bool in_graph = false) const {
    pvector<SGOffset> offsets(vp->partition_width+1);
    for (NodeID_ n = vp->start; n < vp->end+1; n++)
      if (in_graph)
        offsets[n-vp->start] = in_index_[n-vp->start] - in_index_[0];
      else
        offsets[n-vp->start] = out_index_[n-vp->start] - out_index_[0];
    return offsets;
  }

  Range<NodeID_> vertices() const {
    return Range<NodeID_>(num_nodes());
  }

 public:
  bool directed_;
  int64_t num_nodes_;
  int64_t num_edges_;
  DestID_** out_index_;
  DestID_*  out_neighbors_;
  DestID_** in_index_;
  DestID_*  in_neighbors_;
  DestID_* buff_neighbor_;
};

#endif  // GRAPH_H_

