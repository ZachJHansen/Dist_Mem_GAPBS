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
  template <typename T_>  struct OldePartition {
    int64_t N;
    int pe, npes;
    T_ start, end, partition_width, max_width;
  

    OldePartition(int64_t num_nodes) : N(num_nodes){
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
      if (partition_width == 0)
        return npes-1;
      int receiver = node / partition_width;
      if (receiver >= npes)
        receiver = npes - 1;
      return receiver;
    }
  

    // Given a node, determine the local position on assigned PE
    T_ local_pos(T_ node) {
      if (partition_width == 0)
        return node;
      int rec = node / partition_width;
      if (rec >= npes)
        return(node - (npes-1)*partition_width);
      else
        return(node % partition_width);
      //return(node - start);
    }
  };
  

  // bounds for partitioning nodes ~evenly across PEs like a naif (final PE gets remainder)
  // Accessing node v on PE p means accessing node (n/k)*p + v in a complete parent array of n nodes and k PEs
  // Similarly, node V in the complete parent array is the V%(n/k) element in the parent array of PE V/(n/k)
  // (Unless pe = npes-1, then V is the V-(n/k)*p element in the npes-1 PE)
  // should this be generically typed??
  // it would be nice if this incorporated the edge case where num nodes < npes, then first numnodes PEs get 1 node
  template <typename T_=int64_t> struct Partition {
    int64_t N; 
    int pe, npes;
    bool small;
    T_ start, end, partition_width, max_width;
  

    Partition() {}
  

    Partition(int64_t num_nodes, bool ignore=false) : N(num_nodes){
      pe = shmem_my_pe();
      npes = shmem_n_pes();
      //  shmem_barrier_all();    why does this cause deadlock?
      if (num_nodes < npes && !ignore) {
        small = true;
        if (pe < num_nodes) {
          partition_width = 1;
          start = pe;
          end = pe+1;
        } else {
          partition_width = 0;
          start = 0;
          end = 0;
        }
        max_width = 1;
      } else {
        small = false;
        partition_width = N/npes;
        max_width = N - (npes-1)*partition_width;
        start = partition_width * pe;
        if (pe == npes-1) {
          end = N;  
        } else {
          end = start + partition_width; 
        }
      }
    }
    
    // Given a node, determine which PE it belongs to
    int recv(T_ node) {
      int receiver;
      if (small) {
        receiver = node;
      } else {
        receiver = node / partition_width;
        if (receiver >= npes) 
          receiver = npes - 1;
      }
      return receiver;
    }
  

    // Given a node, determine the local position on assigned PE
    T_ local_pos(T_ node) {
      if (small) {
        return 0;
      } else {
        int rec = node / partition_width;
        if (rec >= npes)
          return(node - (npes-1)*partition_width);
        else
          return(node % partition_width);
      }
    }
  

    // Given a local position, determine the global node number
    T_ global_pos(T_ local_pos) {
      return(start + local_pos);
    }
    
    // Return the first node assigned to PE p
    T_ first_node(int p) {
      return partition_width*p;
    }
  

    // Return the partition width of PE p
    T_ p_width(int p) {
      if (small) {
        if (pe < N)
          return 1;
        else
          return 0;
      } else {
        if (pe == npes-1)
          return(N - (npes-1)*partition_width);
        else
          return N/npes;
      }
    }
  };
  

  // Note: Previous partition expects zero-indexed NodeIDs
  // So does Round Robin, but it handles them like one-indexed IDs
  template <typename T_> struct RoundRobin {
    int pe, npes;
    T_ *local_width_, *max_width_;
    
    RoundRobin() {
      pe = shmem_my_pe();
      npes = shmem_n_pes();
      local_width_ = (T_*) shmem_malloc(sizeof(T_));
      *local_width_ = 0;
      max_width_ = (T_*) shmem_malloc(sizeof(T_));
      *max_width_ = -1;
    }
  

    // which PE should K be assigned to?
    int owner(T_ k) {
      k++;                                // change to one-indexing
      if (k < 1) {
        printf("Can't partition 1 element! Breaking...\n");
        shmem_global_exit(1);
        exit(1);
      }
      if (k < npes) {
        return(k-1);
      } else {
        if (k % npes == 0)
          return(npes-1);
        else
          return((k % npes)-1);
      }
    }
  

    void inc_local() {
      (*local_width_)++;
    }
  

    T_ local_width() {
      return *local_width_;
    }
  

    T_ max_width() {
      if (*max_width_ == -1) {
        printf("Max width has not been finalized yet!\n");
        return -1;
      } else {
        return *max_width_;      
      }
    }
  

    // This should only be called once partitioning is over
    // Must be called by all PEs
    T_ finalize_max_width() {
      long* temp_pSync = (long*) shmem_calloc(SHMEM_REDUCE_SYNC_SIZE, sizeof(SHMEM_SYNC_VALUE));
      int* temp_pWrk = (int*) shmem_calloc(SHMEM_REDUCE_MIN_WRKDATA_SIZE, sizeof(SHMEM_SYNC_VALUE)); 
      for (int i = 0; i < SHMEM_REDUCE_SYNC_SIZE; i++)
        temp_pSync[i] = SHMEM_SYNC_VALUE;
      for (int i = 0; i < SHMEM_REDUCE_MIN_WRKDATA_SIZE; i++)
        temp_pWrk[i] = SHMEM_SYNC_VALUE;
      shmem_barrier_all();
      if (sizeof(T_) == sizeof(int)) {
      //  printf("Huzzah- max_width_ = %p, local = %d\n", (void*) max_width_, *local_width_);
        shmem_int_max_to_all(max_width_, local_width_, 1, 0, 0, npes, temp_pWrk, temp_pSync);
      } else {
        printf("Haven't finished implementing generically typed RoundRobin\n");
        shmem_global_exit(0);
        exit(0);
      } 
      shmem_barrier_all();
      shmem_free(temp_pSync);
      shmem_free(temp_pWrk);
      return *max_width_;
    }
  };
   
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
  

/*    void push_back(T_ val) {
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
  */

  // If PEs were to ever call push_back on a symmetric pvector, they would have to call this then push_back_fence
  // But I think this is unneccessary, because while reading edge lists PEs push_back to local pvectors
  void push_back(T_ val) {
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


