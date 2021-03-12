// Copyright (c) 2015, The Regents of the University of California (Regents)
// See LICENSE.txt for license details

#ifndef GENERATOR_H_
#define GENERATOR_H_

#include <cmath>
#include <algorithm>
#include <cinttypes>
#include <iostream>
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
    el.set_combined_length(num_edges_);
    Partition<int64_t> block_partition(num_blocks, true);
    shmem_barrier_all();
    for (int64_t block = (block_size*block_partition.start); block < (block_size*block_partition.end); block += block_size) {
      rng.seed(kRandSeed + block/block_size);
      for (int64_t e=block; e < std::min(block+block_size, num_edges_); e++) {
        el[e - (block_partition.start * block_size)] = Edge(udist(rng), udist(rng));
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
    EdgeList el(edge_count);                            // partitioned edge list
    el.set_combined_length(num_edges_);                 // complete size of edge list
    Partition<int64_t> block_partition(num_blocks, true);
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
  
  static void InsertWeightsEdgeList(pvector<WEdge> &el) {
    std::mt19937 rng;
    std::uniform_int_distribution<int> udist(1, 255);
    int64_t num_edges = el.combined_length();
    int64_t num_blocks = num_edges / block_size;
    if (num_edges % block_size != 0)
      num_blocks++;
    //printf("num_edges %d block size %d num_blocks %lu\n", num_edges, block_size, num_blocks);
    Partition<int64_t> bp(num_blocks);
    int* weights = (int *) shmem_calloc(bp.max_width * block_size, sizeof(int));                        // Allocate, on each PE, enough space for ~num_blocks/npes * block size weights
    //printf("PE %d | num blocks = %d, bp.max_width = %d, bp.part_width = %d, block_size = %d", bp.pe, num_blocks, bp.max_width, bp.partition_width, block_size);
    for (int64_t block=bp.start; block < bp.max_width+bp.start; block++) {                           // Divide processing of blocks naively between PEs
      //printf("PE %d | start = %d\n, end = %d\n", bp.pe, bp.start, bp.max_width+bp.start);
      rng.seed(kRandSeed + block);
      //printf("PE %d | SEEd: %d\n", bp.pe, kRandSeed+block);
      int64_t local_block = (block - bp.start) * block_size;                                              // Calculate local block ID (block is the global block ID within a non-partitioned array)
      //printf("PE %d block %d local block %d\n", bp.pe, block, local_block);
      for (int64_t e=local_block; e < local_block+block_size; e++) {                          // Fill symmetric weights array using local block ID, seeding based on global block ID
        weights[e] = udist(rng);
      }
    }
    /*printf("PE %d ", bp.pe);
    for (int i = 0; i < bp.max_width*block_size; i++)
      printf("%d ", weights[i]);
    printf("\n");*/
    Partition<int64_t> ep(num_edges);                                      // 
    for (int64_t i = 0; i < el.size(); i++) {                                           // Update local region of edge list from weights stored in symmetric memory
      int64_t global_edge = ep.npes * i + ep.pe;                                                            // Position of round-robin partitioned edge i within complete edge list
      int owner = -1;      
      int64_t count = 0;
      for (int k = 0; k < bp.npes-1; k++) {                                       // Find the first PE where the weight associated with edge global_edge is saved
        if (count <= global_edge && global_edge < (count + bp.partition_width*block_size)) {
          owner = k;
          break;
        }
        //count += bp.max_width * block_size;
        count += bp.partition_width * block_size;
      }
      if (owner == -1) 
        owner = bp.npes-1;
      //printf("PE %d | edge %lu has global pos %lu with local pos %lu and owner %d\n", bp.pe, i, global_edge, global_edge-count, owner);
      //shmem_get64(dest, weights+(global_edge-count), 1, owner);
      el[i].v.w = static_cast<WeightT_>(shmem_int_g(weights+(global_edge-count), owner));
      //el[i].v.w = *dest;
    }
    /*for (auto e : el)
      printf("PE %d) = %d ", bp.pe, e.v.w);
    printf("\n");*/
  }

  // Overwrites existing weights with random from [1,255]
  static void InsertWeights(pvector<WEdge> &el, int option) {
    if (option == 0) {
      InsertWeightsEdgeList(el);
    } else {
      InsertWeightsSynthetic(el, option);
    }
  }

  static void InsertWeightsSynthetic(pvector<WEdge> &el, int option) {
    int pe = shmem_my_pe();
    int npes = shmem_n_pes();
    int64_t num_edges = el.combined_length();
    std::mt19937 rng;
    std::uniform_int_distribution<int> udist(1, 255);
    // calculate el position within combined el
    NodeID_ el_offset = 0;
    int64_t num_blocks = num_edges / block_size;
    if (num_edges % block_size != 0)
      num_blocks++;
    if (option == 1) { // if el was generated, then pes 0 - (npes-2) got n*block_size edges, PE npes-1 got remainder
      int64_t blocks_per_pe = num_blocks / npes;
      int i = 0;
      while (i < pe){
        el_offset += blocks_per_pe * block_size;
        i++;
      }
    } else {
      printf("Other graph formats not yet supported\n");
      shmem_global_exit(0);
      exit(1);
    }
    // determine the block the edge would belong to
    int64_t previous_block_count = el_offset / block_size;              // How many complete blocks belonged to previous el partitions?
    int64_t block = previous_block_count * block_size;                  // where does this block begin?
    if (false/*el_offset % block_size == 0*/) {                                  // Edge partitions are aligned with block boundaries
      while ((block-el_offset) < el.size()) {
        rng.seed(kRandSeed + block/block_size);
        if (block < el_offset)                                          // Other PEs have processed previous edges belonging to this block
          rng.discard(el_offset - block);
        for (int64_t e = block; e < std::min(block+block_size, (int64_t) el.size()); e++) {
          el[e - el_offset].v.w = static_cast<WeightT_>(udist(rng));
        }      
        block += block_size;
      }
    } else {                                                            // Edge partition splits the block
      rng.seed(kRandSeed + block/block_size);
      if (block < el_offset)                                          // Other PEs have processed previous edges belonging to this block
        rng.discard(el_offset - block);
      int64_t e = el_offset;
      if (el_offset + el.size() > block + block_size) {                                         // Edge list extends past this block's boundary
        int64_t remainder = (block+block_size) - el_offset;
        for (int64_t r = 0; r < remainder; r++) {                                               // Process remaining portion of block with current seed 
          el[e - el_offset].v.w = static_cast<WeightT_>(udist(rng));
          e++;
        }
        block += block_size;                                                    // Now e and block should be aligned
        while ((block-el_offset) < el.size()) {
          rng.seed(kRandSeed + block/block_size);
          if (block < el_offset)                                          // Other PEs have processed previous edges belonging to this block
            rng.discard(el_offset - block);
          for (int64_t e = block; e < std::min(block+block_size, (int64_t) el.size()+el_offset); e++) {
            el[e - el_offset].v.w = static_cast<WeightT_>(udist(rng));
          }      
          block += block_size;
        }
      } else {                                                                  // Block extends past EL boundary
        for (int64_t remainder = 0; remainder < el.size(); remainder++) {       // Process entire edge list with current seed 
          el[e - el_offset].v.w = static_cast<WeightT_>(udist(rng));
          e++;
        }
      }
    }
    std::ofstream orig_out;
    int* PRINTER = (int*) shmem_calloc(1, sizeof(int));
    shmem_barrier_all();
    pe = shmem_my_pe();
    shmem_int_wait_until(PRINTER, SHMEM_CMP_EQ, pe);           // wait until previous PE puts your pe # in PRINTER
    orig_out.open("/home/zach/projects/Dist_Mem_GAPBS/Dist_Mem_GAPBS/weights_output.txt", std::ios::app);
    for(WEdge wait : el)
      orig_out << wait.v.w << std::endl;
    orig_out.close();
    if (pe != shmem_n_pes()-1)
      shmem_int_p(PRINTER, pe+1, pe+1);             // who's next?
  }

 private:
  int scale_;
  int64_t num_nodes_;
  int64_t num_edges_;
  static const int64_t block_size = 1<<18;
  //static const int64_t block_size = 3;
};

#endif  // GENERATOR_H_
