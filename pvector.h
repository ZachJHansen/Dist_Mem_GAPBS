// Copyright (c) 2015, The Regents of the University of California (Regents)
// See LICENSE.txt for license details

#ifndef PVECTOR_H_
#define PVECTOR_H_

#include "shmem.h"
#include <algorithm>


/*
GAP Benchmark Suite
Class:  pvector
Author: Scott Beamer

Vector class with ability to not initialize or do initialize in parallel
 - std::vector (when resizing) will always initialize, and does it serially
 - When pvector is resized, new elements are uninitialized
 - Resizing is not thread-safe

pvectors in symmetric memory are constructed when the symmetric boolean is true
the internal arrays should be symmetric to support gets and puts
*/

// bounds for partitioning nodes ~evenly across PEs like a naif (final PE gets remainder)
// Accessing node v on PE p means accessing node (n/k)*p + v in a complete parent array of n nodes and k PEs
// Similarly, node V in the complete parent array is the V%(n/k) element in the parent array of PE V/(n/k)
// (Unless pe = npes-1, then V is the V-(n/k)*p element in the npes-1 PE)
// should this be generically typed??
template <typename T_=int64_t> struct Partition {
  int64_t N; 
  int pe, npes;
  T_ start, end, partition_width, max_width;

  Partition() {}

  Partition(int64_t num_nodes) : N(num_nodes){
    pe = shmem_my_pe();
    npes = shmem_n_pes();
    partition_width = N/npes;
    max_width = N - (npes-1)*partition_width;
    start = partition_width * pe;
    if (pe == npes-1) {
      end = N;  
    } else {
      end = start + partition_width; 
    }
  }
  
  // Given a node, determine which PE it belongs to
  int recv(T_ node) {
    int receiver = node / partition_width;
    if (receiver >= npes) 
      receiver = npes - 1;
    return receiver;
  }

  // Given a node, determine the local position on assigned PE
  T_ local_pos(T_ node) {
    int rec = node / partition_width;
    if (rec >= npes)
      return(node - (npes-1)*partition_width);
    else
      return(node % partition_width);
    //return(node - start);
  }

  // Given a local position, determine the global node number
  T_ global_pos(T_ local_pos) {
    return(start + local_pos);
  }

//  void PrintStats(long* PRINT_LOCK) {
//    shmem_set_lock(PRINT_LOCK);
};


template <typename T_>
class pvector {
 public:
  typedef T_* iterator;

  pvector() : start_(nullptr), end_size_(nullptr), end_capacity_(nullptr) {}

  explicit pvector(size_t num_elements, bool symmetric = false) : symmetric_(symmetric) {
    if (symmetric_) {
      start_ = (T_ *) shmem_calloc(num_elements, sizeof(T_));
      end_size_ = start_ + num_elements;
      end_capacity_ = end_size_;
    } else {
      start_ = new T_[num_elements];
      end_size_ = start_ + num_elements;
      end_capacity_ = end_size_;
    }
  }

 pvector(size_t num_elements, T_ init_val, bool symmetric = false) : pvector(num_elements, symmetric) {
    fill(init_val);
  }

  pvector(iterator copy_begin, iterator copy_end) : pvector(copy_end - copy_begin) {
    #pragma omp parallel for
    for (size_t i=0; i < capacity(); i++)
      start_[i] = copy_begin[i];
  }

  // don't want this to be copied, too much data to move
  pvector(const pvector &other) = delete;

  // prefer move because too much data to copy
  pvector(pvector &&other)
      : symmetric_(other.symmetric_), start_(other.start_), end_size_(other.end_size_),
        end_capacity_(other.end_capacity_) {
    other.start_ = nullptr;
    other.end_size_ = nullptr;
    other.end_capacity_ = nullptr;
  }

  // want move assignment
  pvector& operator= (pvector &&other) {
    //printf("move\n");
    symmetric_ = other.symmetric_;
    start_ = other.start_;
    end_size_ = other.end_size_;
    end_capacity_ = other.end_capacity_;
    other.start_ = nullptr;
    other.end_size_ = nullptr;
    other.end_capacity_ = nullptr;
    return *this;
  }

  ~pvector() {
    if (start_ != nullptr) {
      if (symmetric_) {
        shmem_free(start_);
      } else {
        delete[] start_;
      }
    }
  }

  // not thread-safe
  void reserve(size_t num_elements) {
    if (num_elements > capacity()) {
      if (symmetric_) { 
        T_ *new_range = (T_ *) shmem_calloc(num_elements, sizeof(T_));
        #pragma omp parallel for
        for (size_t i=0; i < size(); i++)
          new_range[i] = start_[i];
        end_size_ = new_range + size();
        shmem_free(start_);
        start_ = new_range;
        end_capacity_ = start_ + num_elements;
        shmem_barrier_all();
      } else {
        T_ *new_range = new T_[num_elements];
        #pragma omp parallel for
        for (size_t i=0; i < size(); i++)
          new_range[i] = start_[i];
        end_size_ = new_range + size();
        delete[] start_;
        start_ = new_range;
        end_capacity_ = start_ + num_elements;
      }
    }
  }

  bool empty() {
    return end_size_ == start_;
  }

  void clear() {
    end_size_ = start_;
  }

  void resize(size_t num_elements) {
    reserve(num_elements);
    end_size_ = start_ + num_elements;
  }

  T_& operator[](size_t n) {
    return start_[n];
  }

  const T_& operator[](size_t n) const {
    return start_[n];
  }

  void push_back(T_ val) {
    if (size() == capacity()) {
      size_t new_size = capacity() == 0 ? 1 : capacity() * growth_factor;     // if capacity == 0 newsize = 1 else newsize = cap*gf
      reserve(new_size);
    }
    if (symmetric_) {
      *end_size_ = val;
      end_size_++;
      shmem_barrier_all();                                      // do we want all pes to push back the same element, or broadcast individual elements?
      //shmem_broadcast
    } else {
      *end_size_ = val;
      end_size_++;
    }
  }

  void fill(T_ init_val) {
    #pragma omp parallel for
    for (T_* ptr=start_; ptr < end_size_; ptr++)
      *ptr = init_val;
    if (symmetric_)
      shmem_barrier_all();
  }

  size_t capacity() const {
    return end_capacity_ - start_;
  }

  size_t size() const {
    return end_size_ - start_;
  }

  iterator begin() const {
    return start_;
  }

  iterator end() const {
    return end_size_;
  }

  T_* data() const {
    return start_;
  }

  void swap(pvector &other) {
    std::swap(start_, other.start_);
    std::swap(end_size_, other.end_size_);
    std::swap(end_capacity_, other.end_capacity_);
  }
  
  // Combines the parent pvector from all PEs into a single pvector, which is returned
  pvector<T_> combine(int num_nodes, int pe, int npes, long* pSync) {
    if (!symmetric_) {
      printf("Can't combine pvectors that don't occur in symmetric memory!\n");
      shmem_global_exit(1);
      exit(1);
    } else {
      int start, end;
      int offset = num_nodes/npes;
      start = offset * pe;
      if (pe == npes-1) {
        end = num_nodes;  
      } else {
        end = start + offset; 
      }
      T_* src = (T_ *) shmem_calloc(offset + npes - 1, sizeof(T_));        // Max number of elems any pe can have (offset + max remainder)
      pvector<T_> dest(num_nodes, true);
      #pragma omp parallel for
      for (int n = start; n < end; n++) {
        if (start_[n-start] < -1) {
          src[n-start] = -1;
        } else {
          src[n-start] = start_[n-start];
        }
      }
      if (sizeof(T_) <= 32) {
        shmem_collect32(dest.begin(), src, end-start, 0, 0, npes, pSync);
      } else if (sizeof(T_) <= 64) {                                                                     // else pray it fits in 64 bits?
        shmem_collect64(dest.begin(), src, end-start, 0, 0, npes, pSync);
      } else {
        printf("Requested type for NodeID is larger than 64 bits! pvector.h -> combine method cannot be used. Giving up.\n");
        shmem_global_exit(1);
        exit(1);
      }
      return dest;
    }
  }

  // currently broken. do we even want to recombine? isnt the point of partitioning the node list is too long?
  /*pvector<T_> combine(Partition<T_> vp, long* pSync, bool bfs = false) {
    if (!symmetric_) {
      printf("Can't combine pvectors that don't occur in symmetric memory!\n");
      shmem_global_exit(1);
      exit(1);
    } else {
      pvector<T_> dest(vp.N, true);
      if (bfs) {
        for (T_ n = vp.start; n < vp.end; n++) {
          if (start_[vp.local_pos(n)] < -1) 
            start_[vp.local_pos(n)] = -1;
        }
      }
      printf("Address: %p\n", (void *) start_);
      if (sizeof(T_) <= 32) {
        shmem_collect32(dest.begin(), start_, vp.end-vp.start, 0, 0, vp.npes, pSync);
      } else if (sizeof(T_) <= 64) {                                                                     // else pray it fits in 64 bits?
        shmem_collect64(dest.begin(), start_, vp.end-vp.start, 0, 0, vp.npes, pSync);
      } else {
        printf("Requested type for NodeID is larger than 64 bits! pvector.h -> combine method cannot be used. Giving up.\n");
        shmem_global_exit(1);
        exit(1);
      }
    }
  }*/

 private:
  T_* start_;
  T_* end_size_;
  T_* end_capacity_;
  bool symmetric_;
  static const size_t growth_factor = 2;
};

#endif  // PVECTOR_H_
