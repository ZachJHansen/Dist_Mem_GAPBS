  // Copyright (c) 2015, The Regents of the University of California (Regents)
  // See LICENSE.txt for license details

  #ifndef PARTITION_H_
  #define PARTITION_H_
  
  #include "shmem.h"
  
  /*
  DMM-GAPBS
  Class:  partition
  Author: Zach Hansen
  - Structs for partitioning the graph across PEs
  - Currently supports naive and round robin partitioning 
    (with slightly different alternatives for edge cases such as N < npes) 
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
      //  shmem_barrier_all();    why does this cause deadlock? b/c PEs shouldn't need everyone to calculate these bounds, so they call it alone
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
   
  #endif  // PARTITION_H_



