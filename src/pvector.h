  // Copyright (c) 2015, The Regents of the University of California (Regents)
  // See LICENSE.txt for license details
  
  #ifndef PVECTOR_H_
  #define PVECTOR_H_
  
  #include "partition.h"
  #include <algorithm>
  
  /*
  GAP Benchmark Suite
  Class:  pvector
  Author: Scott Beamer
  
  Vector class with ability to not initialize or do initialize in parallel
   - std::vector (when resizing) will always initialize, and does it serially
   - When pvector is resized, new elements are uninitialized
   - Resizing is not thread-safe
  */

  /*
  DMM-GAPBS
  Author: Zach Hansen  
  Adaptation Notes:
   - Supports traditional and symmetric heap pvectors
   - pvectors in symmetric memory are constructed when the symmetric boolean is true,
     the internal arrays should be symmetric to support gets and puts
   - Typically symmetric pvector contents are not synchronized, and do not need to resize
  */
   
  template <typename T_>
  class pvector {
   public:
    typedef T_* iterator;
  

    pvector() : start_(nullptr), end_size_(nullptr), end_capacity_(nullptr) {}
  

    explicit pvector(size_t num_elements, bool symmetric = false) : symmetric_(symmetric) {
      pe = shmem_my_pe();
      npes = shmem_n_pes();
      local_width_ = (long*) shmem_calloc(1, sizeof(long));
      resize_flag_ = (long *) shmem_calloc(1, sizeof(long));  // 0 means no PE needs to resize, 1 means a PE is requesting a resize
      fence_counter_ = (int *) shmem_calloc(1, sizeof(int));      // the number of PEs that have escaped a pushback region
      if (symmetric_) {
        if (num_elements == 0)
          start_ = (T_ *) shmem_calloc(1, sizeof(T_));                            // is this going to cause problems? callocing with 0 elems returns nullptr
        else 
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
      combined_size_ = other.combined_size_;
      local_width_ = other.local_width_;
      max_width_ = other.max_width_;
      symmetric_ = other.symmetric_;
      start_ = other.start_;
      end_size_ = other.end_size_;
      end_capacity_ = other.end_capacity_;
      other.start_ = nullptr;
      other.end_size_ = nullptr;
      other.end_capacity_ = nullptr;
      other.local_width_ = nullptr;
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
          if (!new_range) {
            printf("PE %d received nullptr from pvector allocation in reserve with %lu neighbors\n", pe, num_elements);
            shmem_global_exit(1);
            exit(1);
          }
          #pragma omp parallel for
          for (size_t i=0; i < size(); i++)
            new_range[i] = start_[i];
          end_size_ = new_range + size();
          shmem_free(start_);
          start_ = new_range;
          end_capacity_ = start_ + num_elements;
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
  

    bool symmetric() const {
      return symmetric_;
    }
  

    void clear() {
      end_size_ = start_;
    }
  

    void resize(size_t num_elements) {
      reserve(num_elements);
      end_size_ = start_ + num_elements;
    }
  

    // Each PE only holds a portion of an edge list typically
    // So it's useful to remember the combined el size to avoid broadcasts 
    void set_combined_length(int64_t len) {
      combined_size_ = len;
    }
  

    int64_t combined_length() {
      return combined_size_;
    }
  

    // i know its a narrowing conversion but i really need v1 to work...
    void set_widths(size_t max_width, size_t local_width) {
      *local_width_ = (long) local_width;
      max_width_ = max_width;
      end_size_ = start_ + local_width;
    }
  

    T_ max_width() {
      return max_width_;
    }
  
    // Get the local partition width of PE k
    long local_width(int k) const {
      if (k == pe) {
        return *local_width_;
      } else {
        long lw = shmem_long_g(local_width_, k);
        return lw;  
      }
    }


    T_& operator[](size_t n) {
      return start_[n];
    }
  

    const T_& operator[](size_t n) const {
      return start_[n];
    }
  

    void push_back(T_ val) {
      if (symmetric_) {                          
        if (size() == capacity()) {
          size_t new_size = capacity() == 0 ? 1 : capacity() * growth_factor;     // if capacity == 0 newsize = 1 else newsize = cap*gf
          reserve(new_size);
        }
        *end_size_ = val;
        end_size_++;
        shmem_barrier_all();                                      // all pes push back the same element
      } else {
        if (size() == capacity()) {
          size_t new_size = capacity() == 0 ? 1 : capacity() * growth_factor;     // if capacity == 0 newsize = 1 else newsize = cap*gf
          reserve(new_size);
        }
        *end_size_ = val;
        end_size_++;
      }
    }

  // If PEs were to ever call push_back on a symmetric pvector, they would have to call this then push_back_fence
  // But I think this is unneccessary, because while reading edge lists PEs push_back to local pvectors
  void symmetric_push_back(T_ val) {
    if (symmetric_) {                          
      if (shmem_test_lock(resize_flag_) == 1) {                                 // a foreign PE has run out of space and is requesting a resize
        size_t new_size = capacity() == 0 ? 1 : capacity() * growth_factor;    
        reserve(new_size);
        *end_size_ = val;
        end_size_++;
      } else {
        if (size() == capacity()) {                                             // the calling PE needs to request a resize
          shmem_set_lock(resize_flag_);                                         // alerts other PEs that a resize is required
          if (size() < capacity()) {                                            // corner case: two PEs request a resize simultaneously, one acquires lock and tries to redundantly resize
            shmem_clear_lock(resize_flag_);
          } else {
            size_t new_size = capacity() == 0 ? 1 : capacity() * growth_factor; // if capacity == 0 newsize = 1 else newsize = cap*gf
            reserve(new_size);
            shmem_clear_lock(resize_flag_);
          }
        }
        *end_size_ = val;
        end_size_++;
      }
    } else {
      if (size() == capacity()) {
        size_t new_size = capacity() == 0 ? 1 : capacity() * growth_factor;     // if capacity == 0 newsize = 1 else newsize = cap*gf
        reserve(new_size);
      }
      *end_size_ = val;
      end_size_++;
    }
  }


  // This method must be called after a region where PEs push back elements to symmetric pvectors
  // If some PEs are done pushing back, but others are not and need help resizing, this method avoids the deadlock
  void push_back_fence() {
    int arrived = shmem_int_atomic_fetch_inc(fence_counter_, 0) + 1;            // register arrival at fence with counter on PE 0
    while (arrived < npes) {
      if (shmem_test_lock(resize_flag_) == 1) {                                 // a foreign PE has run out of space and is requesting a resize
        size_t new_size = capacity() == 0 ? 1 : capacity() * growth_factor;    
        reserve(new_size);
      }
      arrived = shmem_int_atomic_fetch(fence_counter_, 0);
    }
    shmem_barrier_all();
    if (pe == 0)
      *fence_counter_ = 0;                                                      // reset the counter
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
    pvector<long> combine(int64_t num_nodes, long* pSync) {
      if (!symmetric_) {
        printf("Can't combine pvectors that don't occur in symmetric memory!\n");
        shmem_global_exit(1);
        exit(1);
      } else {
        Partition<int> vp(num_nodes);
        pvector<long> dest(num_nodes, true);
        shmem_collect64(dest.begin(), start_, vp.end-vp.start, 0, 0, vp.npes, pSync); 
        shmem_barrier_all();
        return dest;
      }
    }
  

    pvector<int> combine(int64_t num_nodes, int pe, int npes, long* pSync) {
      if (!symmetric_) {
        printf("Can't combine pvectors that don't occur in symmetric memory!\n");
        shmem_global_exit(1);
        exit(1);
      } else {
        Partition<int> vp(num_nodes);
        T_* src = (T_ *) shmem_calloc(vp.max_width, sizeof(T_));        // Max number of elems any pe can have (offset + max remainder)
        pvector<int> dest(num_nodes, true);
        #pragma omp parallel for
        for (int n = vp.start; n < vp.end; n++) {
          if (start_[n-vp.start] < -1) {
            src[n-vp.start] = -1;
          } else {
            src[n-vp.start] = start_[n-vp.start];
          }
        }
        shmem_collect32(dest.begin(), src, vp.end-vp.start, 0, 0, npes, pSync);
        shmem_barrier_all();
        return dest;
      }
    }
  

   private:
    T_* start_;
    T_* end_size_;
    T_* end_capacity_;
    int64_t combined_size_;
    long* local_width_;
    size_t max_width_;
    bool symmetric_;
    long* resize_flag_;
    int* fence_counter_;
    int pe;
    int npes;
    static const size_t growth_factor = 2;
  };
  

  #endif  // PVECTOR_H_



