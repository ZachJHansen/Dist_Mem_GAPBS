// Copyright (c) 2015, The Regents of the University of California (Regents)
// See LICENSE.txt for license details

#ifndef BITMAP_H_
#define BITMAP_H_

#include <algorithm>
#include <cinttypes>

#include "platform_atomics.h"


/*
GAP Benchmark Suite
Class:  Bitmap
Author: Scott Beamer

Parallel bitmap that is thread-safe
 - Can set bits in parallel (set_bit_atomic) unlike std::vector<bool>

New version: build bitmap in symmetric memory
Changes to local bitmaps aren't broadcast to other PEs until merge is called
This assumes that the bitmaps don't ever have their bits unset (1 to 0)
*/


class Bitmap {
 public:
  explicit Bitmap(size_t size, bool symmetric = false) : symmetric_(symmetric) {
    num_words = (size + kBitsPerWord - 1) / kBitsPerWord;
    if (symmetric_) {
      start_ = (uint64_t *) shmem_calloc(num_words, sizeof(uint64_t));
    } else {
      start_ = new uint64_t[num_words]; 
    }
    end_ = start_ + num_words;
  }

  ~Bitmap() {
    if (symmetric_) {
      shmem_free(start_);
    } else {
      delete[] start_; 
    }
  }

  void reset() {
    std::fill(start_, end_, 0);
    shmem_barrier_all();
  }

  // Bitmap of complete adjacency list is not partitioned
  // Graph itself is partitioned, however
  // So each PE maintains a record of which region of the bitmap it's assigned nodes belong to
  // If PE 0 handles nodes 0-2, PE 0 handles Adj[0], Adj[1], Adj[2]
  void init_bitmap_offsets(const Graph &g, Partition<NodeID> vp) {
    NodeID node_interval_start, node_interval_end;
    adj_start_ = (size_t *) shmem_calloc(1, sizeof(size_t));
    size_t ptrdiff = 0;                                 
    for (int i = 0; i < vp.pe; i++) {
      node_interval_start = vp.partition_width * i;
      node_interval_end = vp.partition_width * (i+1) - 1;
      ptrdiff += (g.out_neigh(node_interval_end).finish() - g.out_neigh(node_interval_start).start());
    }
    *adj_start_ = ptrdiff;                                      // Cumulative ptrdiff of adjacency lists on all preceding PEs
    shmem_barrier_all();
  }

  void set_bit(size_t pos) {
    start_[word_offset(pos)] |= ((uint64_t) 1l << bit_offset(pos));
  }

  void set_bit_atomic(size_t pos) {
    uint64_t old_val, new_val;
    do {
      old_val = start_[word_offset(pos)];
      new_val = old_val | ((uint64_t) 1l << bit_offset(pos));
    } while (!compare_and_swap(start_[word_offset(pos)], old_val, new_val));
  }

  // Exists for one use case: v is in u's adjacency list, but the
  // complete adjacency list is partitioned across PEs. The bitmap is
  // not partitioned, so <u,v>'s location in a complete adj list must be reconstructed
  // Each PE maintains its own bitmap, so setting a bit is necessarily atomic (no PE accesses someone elses bitmap)
  // Eventually they are merged with an or to all
  void set_bit_partitioned(const Graph &g, NodeID u, NodeID &v, Partition<NodeID> vp) {
    size_t u_start = shmem_size_g(adj_start_, vp.recv(u));
    size_t local_size = &v - g.out_neigh(vp.first_node(vp.recv(u))).start();
    size_t pos = u_start + local_size;
    start_[word_offset(pos)] |= ((uint64_t) 1l << bit_offset(pos));
  }

  bool get_bit(size_t pos) const {
    return (start_[word_offset(pos)] >> bit_offset(pos)) & 1l;
  }

  // so much communication overhead for one bit...
  bool get_bit_partitioned(const Graph &g, NodeID u, NodeID &v, Partition<NodeID> vp) {
    size_t u_start = shmem_size_g(adj_start_, vp.recv(u));
    size_t local_size = &v - g.out_neigh(vp.first_node(vp.recv(u))).start();
    size_t pos = u_start + local_size;
    return (start_[word_offset(pos)] >> bit_offset(pos)) & 1l;
  }

  void swap(Bitmap &other) {
    std::swap(start_, other.start_);
    std::swap(end_, other.end_);
  }

  // Each PE has a local copy of the bitmap, the set bits on each bitmap need to be combined with bitwise OR
  void merge(long long *pwrk, long *pSync)  {
    // Is a longlong going to work the same as a uinst64_t?
    // I think yes: the bits are all we care about, not the decimal value
    /// So signed vs unsigned shouldnt matter, as long as each word is 64 bits long
    if (symmetric_) {
      shmem_longlong_or_to_all((long long *) start_, (long long *) start_, num_words, 0, 0, shmem_n_pes(), pwrk, pSync);
    } else {
      printf("Bitmaps that do not exist in symmetric memory do not support merging\n");
    }
  }

 private:
  bool symmetric_;
  uint64_t *start_;
  uint64_t *end_;
  size_t *adj_start_;
  uint64_t num_words;
  static const uint64_t kBitsPerWord = 64;
  static uint64_t word_offset(size_t n) { return n / kBitsPerWord; }
  static uint64_t bit_offset(size_t n) { return n & (kBitsPerWord - 1); }
};

#endif  // BITMAP_H_
