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
*/

/*
DMM-GAPBS
Author: Zach Hansen
Adaptation Nodes:
 - Partitioned Sliding Queue (PEs replace threads)
 - Every PE maintains a QueueBuffer. When flushing, contents are distributed round robin to all PEs
 - If the requested PE does not have space, then the next PE with space is given the value 
 - Window on a PE is incremented when an element is distributed to it
 - Sliding Queue is symmetric but unsynched. 
 - At any time, the concatenation of all visible PE windows should be the same set of values as a non-partitioned window
 - During kernel execution, each PE processes the contents of their Queue
*/

template <typename T>
class QueueBuffer;

template <typename T>
class SlidingQueue {
  T* shared;
  size_t* shared_in;
  size_t* shared_out_start;
  size_t* shared_out_end;
  size_t shared_size_;
  long* QLOCK_;
  friend class QueueBuffer<T>;

 public:
  explicit SlidingQueue(size_t shared_size, long* QLOCK) : QLOCK_(QLOCK), shared_size_(shared_size) {
    shared = (T *) shmem_calloc(shared_size, sizeof(T));
    if (!shared) {
      printf("Allocating a sliding queue of size %lu failed. Try increasing symmetric heap size.\n", shared_size);
      shmem_global_exit(1);
      exit(1); 
    }
    shared_in = (size_t *) shmem_malloc(sizeof(size_t));
    shared_out_start = (size_t *) shmem_malloc(sizeof(size_t));
    shared_out_end = (size_t *) shmem_malloc(sizeof(size_t));
    reset();
  }

  ~SlidingQueue() {
    shmem_free(shared);
    shmem_free(shared_in);
    shmem_free(shared_out_start);
    shmem_free(shared_out_end);
  }

  // Push an element back directly to the calling PE's queue
  void push_back(T to_add) {
    shared[(*shared_in)++] = to_add;                                                   
  }

  // Search for a PE with a non-empty queue
  bool empty() const {
    bool flag = true;
    for (int i = 0; i < shmem_n_pes(); i++) {
      if (shmem_size_g(shared_out_start, i) != shmem_size_g(shared_out_end, i)) {
        flag = false;
        break;
      }
    }
    return flag;
  }

  void reset() {                                // Sync point
    shmem_barrier_all();
    *shared_out_start = 0;
    *shared_out_end = 0;
    *shared_in = 0;
    shmem_barrier_all();
  }

  // Each PE maintains different window boundaries, so windows are slid locally
  void slide_window() {                         // Sync point
    shmem_barrier_all();
    *shared_out_start = *shared_out_end;
    *shared_out_end = *shared_in;
    shmem_barrier_all();
  }

  typedef T* iterator;

  iterator begin() const {
    return shared + *shared_out_start;
  }

  iterator end() const {
    return shared + *shared_out_end;
  }

  // Size of the unpartitioned window
  size_t size() const {
    size_t size = 0;
    for (int i = 0; i < shmem_n_pes(); i++) 
      size += shmem_size_g(shared_out_end, i) - shmem_size_g(shared_out_start, i);
    return size;
  }
};


template <typename T>
class QueueBuffer {
  size_t in;                    // How many incoming elements have been added to QB
  int* scatter_counters;        // How many elements have been sent to each PE during a flush
  T* local_queue;       
  size_t* foreign_windows;          // Addresses on each PE where window begins
  SlidingQueue<T> &sq;
  const size_t local_size;
  size_t data_size;
  int pe;
  int npes; 

 public:
  explicit QueueBuffer(SlidingQueue<T> &master, size_t given_size = 16384)
                      : sq(master), local_size(given_size) {
    data_size = sizeof(T);
    pe = shmem_my_pe();
    npes = shmem_n_pes();
    in = 0;
    local_queue = new T[local_size];
    scatter_counters = new int[npes];
    foreign_windows = new size_t[npes];
  }

  ~QueueBuffer() {
    delete[] local_queue;
    delete[] scatter_counters;
    delete[] foreign_windows;
  }

  void push_back(T to_add) {
    if (in == local_size)
      flush();
    local_queue[in++] = to_add;
  }


  void flush() {
    if (in != 0) {                                                                              // Avoid accidental flushes of empty buffers
      shmem_set_lock(sq.QLOCK_);                                                                // Lock critical region (the frontier) to avoid simultaneous flushes
      for (int i = 0; i < npes; i++) {                                                          // Get symmetric heap offset where each window begins
        foreign_windows[i] = shmem_size_g(sq.shared_in, i);
        scatter_counters[i] = 0;
      }
      T* dest;
      int current_pe = 0;
      for (size_t i = 0; i < in; i++) {                                                         // Distribute queue buffer contents round robin
        dest = sq.shared + (foreign_windows[current_pe] + scatter_counters[current_pe]); 
        while (sq.shared_size_ <= foreign_windows[current_pe] + scatter_counters[current_pe]) { // While there is not available queue space on the current PE
          current_pe++;                                                                         // search for a PE with space
          if (current_pe == npes) {
            printf("The sliding queue does not have enough space to accomodate all elements.\n");
            shmem_global_exit(1);
            exit(1);
          }
          dest = sq.shared + (foreign_windows[current_pe] + scatter_counters[current_pe]);
        }
        shmem_putmem(dest, local_queue+i, data_size, current_pe);
        scatter_counters[current_pe]++;
        current_pe++;
        if (current_pe >= npes)
          current_pe = 0;
      } 
      for (int i = 0; i < npes; i++) {                                                          // Adjust all PE windows to reflect the number of added elements
        shmem_size_p(sq.shared_in, foreign_windows[i] + scatter_counters[i], i); 
      }
      in = 0;
      shmem_clear_lock(sq.QLOCK_);
    }
  }        
};

#endif  // SLIDING_QUEUE_H_
