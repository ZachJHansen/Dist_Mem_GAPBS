// Copyright (c) 2015, The Regents of the University of California (Regents)
// See LICENSE.txt for license details

#ifndef PVECTOR_H_
#define PVECTOR_H_

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
      size_t new_size = capacity() == 0 ? 1 : capacity() * growth_factor;
      reserve(new_size);
    }
    *end_size_ = val;
    end_size_++;
  }

  void fill(T_ init_val) {
    #pragma omp parallel for
    for (T_* ptr=start_; ptr < end_size_; ptr++)
      *ptr = init_val;
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
  // Assumes nodeids are 32 bit integers
  pvector<int32_t> combine(int num_nodes, int pe, int npes, long* pSync) {
    if (!symmetric_) {
      printf("Can't combine pvectors that don't occur in symmetric memory!\n");
      shmem_global_exit(1);
    } else {
      int start, end;
      int offset = num_nodes/npes;
      start = offset * pe;
      if (pe == npes-1) {
        end = num_nodes;  
      } else {
        end = start + offset; 
      }
      int32_t* src = (int32_t *) shmem_calloc(offset + npes - 1, sizeof(int32_t));        // Max number of elems any pe can have (offset + max remainder)
      pvector<int32_t> dest(num_nodes, true);
      #pragma omp parallel for
      for (int n = start; n < end; n++) {
        if (start_[n-start] < -1) {
          src[n-start] = -1;
        } else {
          src[n-start] = start_[n-start];
        }
      }
      shmem_collect32(dest.begin(), src, end-start, 0, 0, npes, pSync);
      return dest;
    }
  }

 private:
  T_* start_;
  T_* end_size_;
  T_* end_capacity_;
  bool symmetric_;
  static const size_t growth_factor = 2;
};

#endif  // PVECTOR_H_