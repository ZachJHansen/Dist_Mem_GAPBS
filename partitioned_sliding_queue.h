// Copyright (c) 2015, The Regents of the University of California (Regents)
// See LICENSE.txt for license details

#ifndef SLIDING_QUEUE_H_
#define SLIDING_QUEUE_H_

#include <algorithm>
#include "shmem.h"
#include "platform_atomics.h"


/*
GAP Benchmark Suite
Class:  SlidingQueue
Author: Scott Beamer

Double-buffered queue so appends aren't seen until SlideWindow() called
 - Use QueueBuffer when used in parallel to avoid false sharing by doing
   bulk appends from thread-local storage

  Reworked in such a way that an individual PE does the work of a thread
  When the sliding queue in symmetric memory is updated by one PE, other PEs must wait to access
  I really don't want to partition the sliding queue, so to save space it will wrap around
  If it grows past the max size the shared start will be st to 0, shared end will be 0 + size
  copy elements from end of q to beginning, overwriting no-longer-used elems
*/


template <typename T>
class QueueBuffer;

template <typename T>
class SlidingQueue {
  int max;
  T *shared;
  bool wrap;
  size_t shared_in;
  size_t shared_out_start;
  size_t shared_out_end;
  int pe;
  int npes;
  friend class QueueBuffer<T>;

 public:
  explicit SlidingQueue(size_t shared_size) : max(1048576) {
    if (shared_size > max) {
      wrap = true;
      shared = (T *) shmem_calloc(max, sizeof(T));
    } else {
      wrap = false;
      shared = (T *) shmem_calloc(shared_size, sizeof(T));
    }
    reset();
  }

  ~SlidingQueue() {
    shfree(shared);
  }

  // every PE performs the same wrap when pushing back the same element
  void wrapper(T to_add) {
    if (shared_out_start == start) {
      printf("Queue size has exceeded maximum allowance for wrapping! Resizing...\n");
      //shmem_realloc(start, sizeof(T)*(max + 16384);   can't go in crit region
    }
    std::copy(shared_out_start, shared_in, shared);
    shared_in -= shared_out_start;
    shared_out_end -= shared_out_start;
    shared_out_start = start;
    shared[shared_in++] = to_add;
    shmem_clear_lock(WRAP_LOCK);
  }

  // wrapping required during a flush, calling PE
  // wraps everyone's sliding queues
  void flush_wrapper() {
    std::copy(shared_out_start, shared_in, shared);
    *shared_in -= *shared_out_start;
    *shared_out_end -= *shared_out_start;
    *shared_out_start = start;
    shared[(*shared_in)++] = to_add;
    for (int i = 0; i < npes; i++) {
      if (i != pe) {
        shmem_putmem(start, start, sizeof(T)*shared_in, i);
        shmem_size_put(shared_in, shared_in, 1, i);
        shmem_size_put(shared_out_end, shared_out_end, 1, i);
        shmem_size_put(shared_out_start, shared_out_end, 1, i);
      }
    }
  }

  // pushback on a symmetric sliding queue should only occur if all pes are 
  // pushing the same to_add at the same step
  void push_back(T to_add) {
    shmem_barrier_all();
    if (wrap) {
      if (shared_in++ >= max && shmem_test_lock(WRAP_LOCK) == 0) {
        wrapper(to_add);
      } else {
        shared[shared_in++] = to_add;
      }
    } else {
      shared[shared_in++] = to_add;
    }
    shmem_barrier_all();
  }

  bool empty() const {
    return shared_out_start == shared_out_end;
  }

  void reset() {                                // Sync point
    shmem_barrier_all();
    shared_out_start = 0;
    shared_out_end = 0;
    shared_in = 0;
    shmem_barrier_all();
  }

  void slide_window() {                         // Sync point
    shmem_barrier_all();
    shared_out_start = shared_out_end;
    shared_out_end = shared_in;
    shmem_barrier_all();
  }

  typedef T* iterator;

  iterator begin() const {
    return shared + shared_out_start;
  }

  iterator end() const {
    return shared + shared_out_end;
  }

  size_t size() const {
    return end() - begin();
  }
};


template <typename T>
class QueueBuffer {
  size_t in;
  T *local_queue;
  SlidingQueue<T> &sq;
  const size_t local_size;
  size_t data_size;
  int pe;
  int npes; 
  long* QLOCK;

 public:
  explicit QueueBuffer(SlidingQueue<T> &master, long* QL, size_t given_size = 16384)
                      : sq(master), QLOCK(QL), local_size(given_size) {
    data_size = sizeof(T);
    pe = shmem_my_pe();
    npes = shmem_n_pes();
    in = 0;
    local_queue = new T[local_size];
  }

  ~QueueBuffer() {
    delete[] local_queue;
  }

  void push_back(T to_add) {
    if (in == local_size)
      flush();
    local_queue[in++] = to_add;
  }

  void flush() {
    shmem_set_lock(QLOCK);                                                                    // Lock critical region (the frontier) to avoid simultaneous flushes
    if (wrap) {
      if (sq.shared_in - sq.shared_out_start + in) > sq.max/2)) {
        printf("Need to resize queue...\n");
  //      shmem_realloc(sq.shared, sizeof(T)*16384);            cant have this in crit region b/c all pes must participate?
        sq.flush_wrapper();
      } else {
        sq.flush_wrapper();
      }
    }
    if (in != 0) {                                                                              // Avoid accidental flushes of empty buffers
      T *shared_queue = sq.shared;
      size_t copy_start = shmem_ulong_atomic_fetch_add(&(sq.shared_in), in, pe);                // Get start of shared queue incoming region, update local copy of incoming region start position
      size_t copy_end = copy_start + in;                // what if you need to wrap around?
      std::copy(local_queue, local_queue+in, shared_queue+copy_start);                          // Update local copy of shared queue
      for (int i = 0; i < npes; i++){
        if (i != pe){
          shmem_putmem(shared_queue+copy_start, local_queue, data_size*in, i);                  // Update shared queue on all PEs (put a contiguous block of memory determined by data type & size)
          shmem_ulong_put(&(sq.shared_in), &copy_end, (long unsigned) 1, i);                    // Move start of incoming region to end of copied elements
        }
      }
      in = 0;
    }
  }
  shmem_clear_lock(QLOCK);
};

#endif  // SLIDING_QUEUE_H_
