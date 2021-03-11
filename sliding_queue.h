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
    Sliding queues should be symmetric and synched - every PE has the same contents
    Queue buffers are not synched or symmetric
*/


template <typename T>
class QueueBuffer;

template <typename T>
class SlidingQueue {
  T *shared;
  size_t* shared_in;
  size_t shared_out_start;
  size_t shared_out_end;
  long* QLOCK_;
  friend class QueueBuffer<T>;

 public:
  explicit SlidingQueue(size_t shared_size, long* QLOCK) : QLOCK_(QLOCK) {
    shared = (T *) shmem_calloc(shared_size, sizeof(T));
    shared_in = (size_t *) shmem_malloc(sizeof(size_t));
    reset();
  }

  ~SlidingQueue() {
    shmem_free(shared);
    shmem_free(shared_in);
  }

  // Only one PE should call pushback for a given element, 
  // because it updates everyone's copy of the shared sliding queue
  void push_back(T to_add) {
    shmem_set_lock(QLOCK_);
    shared[(*shared_in)++] = to_add;                                                        // update local copy
    for (int i = 0; i < shmem_n_pes(); i++) {                                            // update all PEs with to_add value and incremented shared_in
      if (shmem_my_pe() != i) {
        shmem_int_p(shared+(*shared_in-1), to_add, i);                            // int should be of type T
        shmem_size_p(shared_in, *shared_in, i);
      }
    }
    shmem_clear_lock(QLOCK_);
  }

  bool empty() const {
    return shared_out_start == shared_out_end;
  }

  void reset() {                                // Sync point
    shmem_barrier_all();
    shared_out_start = 0;
    shared_out_end = 0;
    *shared_in = 0;
    shmem_barrier_all();
  }

  void slide_window() {                         // Sync point
    shmem_barrier_all();
    shared_out_start = shared_out_end;
    shared_out_end = *shared_in;
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

 public:
  explicit QueueBuffer(SlidingQueue<T> &master, size_t given_size = 16384)
                      : sq(master), local_size(given_size) {
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
    if (in != 0) {                                                                              // Avoid accidental flushes of empty buffers
      shmem_set_lock(sq.QLOCK_);                                                                    // Lock critical region (the frontier) to avoid simultaneous flushes
      T *shared_queue = sq.shared;
      size_t copy_start = shmem_ulong_atomic_fetch_add(sq.shared_in, in, pe);                // Get start of shared queue incoming region, update local copy of incoming region start position
      size_t copy_end = copy_start + in;
      std::copy(local_queue, local_queue+in, shared_queue+copy_start);                          // Update local copy of shared queue
      for (int i = 0; i < npes; i++){
        if (i != pe){
          shmem_putmem(shared_queue+copy_start, local_queue, data_size*in, i);                  // Update shared queue on all PEs (put a contiguous block of memory determined by data type & size)
          shmem_size_put(sq.shared_in, &copy_end, (size_t) 1, i);                    // Move start of incoming region to end of copied elements
        }
      }
      in = 0;
      shmem_clear_lock(sq.QLOCK_);
    }
  }
};

#endif  // SLIDING_QUEUE_H_
