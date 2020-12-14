// Copyright (c) 2015, The Regents of the University of California (Regents)
// See LICENSE.txt for license details

#ifndef GENERATOR_H_
#define GENERATOR_H_

#include <algorithm>
#include <cinttypes>
#include <random>

#include "graph.h"
#include "pvector.h"
#include "util.h"


/*
GAP Benchmark Suite
Class:  Generator
Author: Scott Beamer

Given scale and degree, generates edgelist for synthetic graph
 - Intended to be called from Builder
 - GenerateEL(uniform) generates and returns the edgelist
 - Can generate uniform random (uniform=true) or R-MAT graph according
   to Graph500 parameters (uniform=false)
 - Can also randomize weights within a weighted edgelist (InsertWeights)
 - Blocking/reseeding is for parallelism with deterministic output edgelist
*/


template <typename NodeID_, typename DestID_ = NodeID_,
          typename WeightT_ = NodeID_>
class Generator {
  typedef EdgePair<NodeID_, DestID_> Edge;
  typedef EdgePair<NodeID_, NodeWeight<NodeID_, WeightT_>> WEdge;
  typedef pvector<Edge> EdgeList;

 public:
  Generator(int scale, int degree) {
    scale_ = scale;
    num_nodes_ = 1l << scale;
    num_edges_ = num_nodes_ * degree;
    if (num_nodes_ > std::numeric_limits<NodeID_>::max()) {
      std::cout << "NodeID type (max: " << std::numeric_limits<NodeID_>::max();
      std::cout << ") too small to hold " << num_nodes_ << std::endl;
      std::cout << "Recommend changing NodeID (typedef'd in src/benchmark.h)";
      std::cout << " to a wider type and recompiling" << std::endl;
      std::exit(-31);
    }
  }

  void PermuteIDs(EdgeList &el) {
    NodeID_ global_e, owner, offset;
    pvector<NodeID_> permutation(num_nodes_);
    std::mt19937 rng(kRandSeed);
    for (NodeID_ n=0; n < num_nodes_; n++)
      permutation[n] = n;
    shuffle(permutation.begin(), permutation.end(), rng);
    for (int64_t e = 0; e < el.size(); e++) {
      el[e] = Edge(permutation[el[e].u], permutation[el[e].v]);
    }
  }

  // Return the number of elements that would be allocated to PE
  // p with NPES n using round robin partitioning
  NodeID_ round_robin_size(NodeID_ num_edges, int p, int n) {
    NodeID_ k = num_edges - (num_edges % n);
    NodeID_ r = num_edges - (n * k);
    int counts[n];
    for (int i = 0; i < n; i++)
      counts[i] = k;
    int j = 0;
    while (j < r) {
      counts[j]++;
      j++;
    }
    return counts[p];
  }

  // All PEs except the last one get X blocks of edges
  // The last PE gets X*block size + the remaining edges
  // So PEs have unbalanced edge lists, reducing the block size should help
  // Splitting blocks across PEs just seems horrible, and the block is used to seed the rng, so...
  EdgeList MakeUniformEL() {
    int pe = shmem_my_pe();
    int npes = shmem_n_pes();
    std::mt19937 rng;
    std::uniform_int_distribution<NodeID_> udist(0, num_nodes_-1);
    int64_t num_blocks = num_edges_ / block_size;
    if (num_edges_ % block_size != 0)
      num_blocks++;
    NodeID_ edge_count;
    int64_t blocks_per_pe = num_blocks / npes;
    if (pe == npes - 1) {
      edge_count = num_edges_ - (blocks_per_pe * (npes - 1)) * block_size;
    } else {
      edge_count = blocks_per_pe * block_size;
    }
    EdgeList el(edge_count);
  //  printf("PE %d | EL creation with %lu edges\n", pe, edge_count);
    Partition<int64_t> block_partition(num_blocks);
    shmem_barrier_all();
    printf("PE %d has bp start %lu, bp end %lu, block size %lu for num_edges_ %lu and num blocks %lu\n", pe, block_partition.start, block_partition.end, block_size, num_edges_, num_blocks); 
    for (int64_t block = (block_size*block_partition.start); block < (block_size*block_partition.end); block += block_size) {
      rng.seed(kRandSeed + block/block_size);
      for (int64_t e=block; e < std::min(block+block_size, num_edges_); e++) {
        //el[0] = Edge(udist(rng), udist(rng));
       // if (e - block_partition.start * block_size < 0) {
    //      ;
         // printf("negative\n");
       // } else {
        el[e - (block_partition.start * block_size)] = Edge(udist(rng), udist(rng));
        //}
      } 
    }
    return el;
  }

  EdgeList MakeRMatEL() {
    int pe = shmem_my_pe();
    int npes = shmem_n_pes();
    const float A = 0.57f, B = 0.19f, C = 0.19f;
    std::mt19937 rng;
    std::uniform_real_distribution<float> udist(0, 1.0f);
    int64_t num_blocks = num_edges_ / block_size;
    if (num_edges_ % block_size != 0)
      num_blocks++;
    NodeID_ edge_count;
    int64_t blocks_per_pe = num_blocks / npes;
    if (pe == npes - 1) {
      edge_count = num_edges_ - (blocks_per_pe * (npes - 1)) * block_size;
    } else {
      edge_count = blocks_per_pe * block_size;
    }
    EdgeList el(edge_count);
    Partition<int64_t> block_partition(num_blocks);
    for (int64_t block = (block_size*block_partition.start); block < (block_size*block_partition.end); block += block_size) {
      rng.seed(kRandSeed + block/block_size);
      for (int64_t e=block; e < std::min(block+block_size, num_edges_); e++) {
        NodeID_ src = 0, dst = 0;
        for (int depth=0; depth < scale_; depth++) {
          float rand_point = udist(rng);
          src = src << 1;
          dst = dst << 1;
          if (rand_point < A+B) {
            if (rand_point > A)
              dst++;
          } else {
            src++;
            if (rand_point > A+B+C)
              dst++;
          }
        }
        el[e - (block_partition.start * block_size)] = Edge(src, dst);
      }
    }
    PermuteIDs(el);
    return el;
  }

  EdgeList GenerateEL(bool uniform) {
    EdgeList el;
    Timer t;
    t.Start();
    if (uniform)
      el = MakeUniformEL();
    else
      el = MakeRMatEL();
    t.Stop();
    PrintTime("Generate Time", t.Seconds());
    return el;
  }

  static void InsertWeights(pvector<EdgePair<NodeID_, NodeID_>> &el, int option) {}

  // Overwrites existing weights with random from [1,255]
  static void InsertWeights(pvector<WEdge> &el, int option) {
    int pe = shmem_my_pe();
    int npes = shmem_n_pes();
    int64_t num_edges = el.combined_size();
    std::mt19937 rng;
    std::uniform_int_distribution<int> udist(1, 255);
    // calculate el position within combined el
    NodeID_ el_offset = 0;
    int64_t num_blocks = num_edges / block_size;
    if (num_edges % block_size != 0)
      num_blocks++;
    // if el was read in from a file, the el was partitioned round robin
    if (option == 0) {
      int i = 0;
      while (i < pe)
        el_offset += round_robin_size(num_edges, i, npes);
    } else if (option == 1) { // if el was generated, then pes 0 - (npes-2) got n*block_size edges, PE npes-1 got remainder
      int64_t blocks_per_pe = num_blocks / npes;
      int i = 0;
      while (i < pe)
        el_offset += blocks_per_pe * block_size;
    } else {
      printf("Other graph formats not yet supported\n");
      shmem_global_exit(0);
      exit(1);
    }
    // determine the block the edge would belong to
    int64_t previous_block_count = el_offset / block_size;              // How many complete blocks belonged to previous el partitions?
    int64_t block = previous_block_count * block_size;                  // where does this block begin?
    bool start = true;
    if (el_offset % block_size == 0) {                                  // Edge partitions are aligned with block boundaries
      while ((block-el_offset) < el.size()) {
        rng.seed(kRandSeed + block/block_size);
        for (int64_t e = block; e < std::min(block+block_size, el.size()); e++) {
          el[e - el_offset].v.w = static_cast<WeightT_>(udist(rng));
        }      
        block += block_size;
      }
    } else {                                                            // Edge partition splits the block
      rng.seed(kRandSeed + block/block_size);
      int64_t e = el_offset;
      for (int64_t remainder = 0; remainder < (el_offset % block_size); remainder++) {  // Process remaining portion of block with initial rng
        el[e - el_offset].v.w = static_cast<WeightT_>(udist(rng));
        e++;
      }
      block += block_size;                                              // Now e and block should be aligned
      while ((block-el_offset) < el.size()) {
        rng.seed(kRandSeed + block/block_size);
        for (int64_t e = block; e < std::min(block+block_size, el.size()); e++) {
          el[e - el_offset].v.w = static_cast<WeightT_>(udist(rng));
        }      
      }
    }
  }

  // Overwrites existing weights with random from [1,255]
/*  static void InsertWeights(pvector<WEdge> &el) {
    #pragma omp parallel
    {
      std::mt19937 rng;
      std::uniform_int_distribution<int> udist(1, 255);
      int64_t el_size = el.size();
      #pragma omp for
      for (int64_t block=0; block < el_size; block+=block_size) {
        rng.seed(kRandSeed + block/block_size);
        for (int64_t e=block; e < std::min(block+block_size, el_size); e++) {
          el[e].v.w = static_cast<WeightT_>(udist(rng));
        }
      }
    }
  }
*/

 private:
  int scale_;
  int64_t num_nodes_;
  int64_t num_edges_;
  static const int64_t block_size = 1<<18;
};

#endif  // GENERATOR_H_
