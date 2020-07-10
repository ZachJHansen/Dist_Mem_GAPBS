// Copyright (c) 2015, The Regents of the University of California (Regents)
// See LICENSE.txt for license details

#include <cinttypes>
#include <limits>
#include <iostream>
#include <queue>
#include <vector>

#include "benchmark.h"
#include "builder.h"
#include "command_line.h"
#include "graph.h"
#include "platform_atomics.h"
#include "pvector.h"
#include "timer.h"


/*
GAP Benchmark Suite
Kernel: Single-source Shortest Paths (SSSP)
Author: Scott Beamer

Returns array of distances for all vertices from given source vertex

This SSSP implementation makes use of the ∆-stepping algorithm [1]. The type
used for weights and distances (WeightT) is typedefined in benchmark.h. The
delta parameter (-d) should be set for each input graph.

The bins of width delta are actually all thread-local and of type std::vector
so they can grow but are otherwise capacity-proportional. Each iteration is
done in two phases separated by barriers. In the first phase, the current
shared bin is processed by all threads. As they find vertices whose distance
they are able to improve, they add them to their thread-local bins. During this
phase, each thread also votes on what the next bin should be (smallest
non-empty bin). In the next phase, each thread copies their selected
thread-local bin into the shared bin.

Once a vertex is added to a bin, it is not removed, even if its distance is
later updated and it now appears in a lower bin. We find ignoring vertices if
their current distance is less than the min distance for the bin to remove
enough redundant work that this is faster than removing the vertex from older
bins.

[1] Ulrich Meyer and Peter Sanders. "δ-stepping: a parallelizable shortest path
    algorithm." Journal of Algorithms, 49(1):114–152, 2003.
*/


// WARNING! Broken C&S

using namespace std;

const WeightT kDistInf = numeric_limits<WeightT>::max()/2;
const size_t kMaxBin = numeric_limits<size_t>::max()/2;
const size_t kBinSizeThreshold = 1000;

inline
void RelaxEdges(const WGraph &g, NodeID u, WeightT delta, pvector<WeightT> &dist, vector <vector<NodeID>> &local_bins, Partition<NodeID> vp) {
  WeightT old_dist, new_dist;
  printf("Relaxing u: %d\n", u);
  //printf("sizeof(WeightT) == %lu\n", sizeof(WeightT));
  //printf("sizoef(int) == %lu\n", sizeof(int));
  for (WNode wn : g.out_neigh(u)) {
    // modifications to local data should also be atomic! 
    shmem_getmem(&old_dist, dist.begin()+vp.local_pos(wn.v), sizeof(WeightT), vp.recv(wn.v));               // recv should be renamed owner, no?
    shmem_getmem(&new_dist, dist.begin()+vp.local_pos(u), sizeof(WeightT), vp.recv(u));               // can a pe fetch data from itself?
    new_dist += wn.w;
    //printf("PE: %d | wn.v = %d | vp.start = %d | vp.end = %d\n", vp.pe, wn.v, vp.start, vp.end);
    printf("PE %d | New dist: %d | old_dist: %d\n", vp.pe, new_dist, old_dist);
    //printf("wn.v = %d | local_pos = %d | recv = %d \n", wn.v, vp.local_pos(wn.v), vp.recv(wn.v));
    while (new_dist < old_dist) {
      /*printf("Dist begins at %p, dist + local_pos(%d) = %p\n", (void *) dist.begin(), wn.v, (void *) (dist.begin()+vp.local_pos(wn.v))); 
      if (shmem_addr_accessible(dist.begin()+vp.local_pos(wn.v), vp.recv(wn.v)))
        printf("C&S at %p on PE %d with cond = %d, val = %d\n", (void *) (dist.begin()+vp.local_pos(wn.v)), vp.pe, old_dist, new_dist);
      */
    //  if (shmem_int_atomic_compare_swap(dist.begin()+vp.local_pos(wn.v), (int) old_dist, (int) new_dist, vp.recv(wn.v))) {              // not sure how to generically type atomics...
      //long x;
      //if ((x = shmem_long_atomic_compare_swap(dist.begin()+vp.local_pos(wn.v), (long) 1, (long) 1, 1)) > 0) {
  //      printf("check! x = %d\n", x);
      if (shmem_int_atomic_fetch(dist.begin()+vp.local_pos(wn.v), vp.recv(wn.v)) == old_dist) {          // these two atomics should really be combined into C&S...
        shmem_int_atomic_swap(dist.begin()+vp.local_pos(wn.v), new_dist, vp.recv(wn.v));
        size_t dest_bin = new_dist/delta;
        if (dest_bin >= local_bins.size())
          local_bins.resize(dest_bin+1);
        local_bins[dest_bin].push_back(wn.v);
        //printf("Check 2\n");
        break;
      }
      if (wn.v >= vp.start && wn.v < vp.end)
        old_dist = dist[wn.v];      // swap failed, recheck dist update & retry
      else
        shmem_getmem(&old_dist, dist.begin()+vp.local_pos(wn.v), sizeof(WeightT), vp.recv(wn.v));               
    }
  }
  printf("No more relaxing\n");
}

pvector<WeightT> Shmem_DeltaStep(const WGraph &g, NodeID source, int delta, long* VOTE_LOCK, long* SINGLE_LOCK) {
  WeightT udist;
  WeightT init = 0;
  WeightT *old_dist, *new_dist;
  Timer t;
  Partition<NodeID> vp(g.num_nodes());
  //printf("PE %d is init dist with width = %lu\n", vp.pe, vp.max_width);
//  pvector<WeightT> dist(vp.max_width, kDistInf, true);                                 // Initialize symmetric partitioned pvector with +inf
  pvector<WeightT> dist(vp.max_width, kDistInf, true);
  //for (int d : dist)
    //printf("Init: %d\n", d);
  if (source >= vp.start && source < vp.end)
    dist[vp.local_pos(source)] = init;
  else
    shmem_putmem(dist.begin()+vp.local_pos(source), &init, sizeof(WeightT), vp.recv(source));
  //printf("PE %d is init frontier with size = %lu\n", vp.pe, g.num_edges_directed());
  pvector<NodeID> frontier(g.num_edges_directed(), true);                               // Symmetric frontier - intentionally overallocated. Do we want a space or performance hit?
  // two element arrays for double buffering curr=iter&1, next=(iter+1)&1
  size_t* shared_indexes = (size_t *) shmem_calloc(2, sizeof(size_t));                  // size_t shared_indexes[2] = {0, kMaxBin};
  shared_indexes[1] = kMaxBin;
  size_t* frontier_tails = (size_t *) shmem_calloc(2, sizeof(size_t));                   // size_t frontier_tails[2] = {1, 0};
  frontier_tails[0] = 1;
  frontier[0] = source;
  t.Start();                                                                    // Timer start and stops are synch points
  vector<vector<NodeID> > local_bins(0);                                      // vector of vector of node ids (thread local?)
  NodeID* iter = (NodeID *) shmem_calloc(1, sizeof(NodeID));                                  // Shared counter for next node (init 0)
  //printf("PE %d | dist.begin = %p, dist.end = %p, front.begin = %p, front.end = %p\n", vp.pe, (void *) dist.begin(),  (void *) dist.end(), (void *) frontier.begin(), (void *) frontier.end());
  while (shared_indexes[(*iter)&1] != kMaxBin) {                                 // if iter is odd, iter bitwise and (&) 1 is 1, if iter is even iter&1 is 0 and enter loop
    //printf("PE %d | While start | iter = %d, curr_front = %lu\n", vp.pe, *iter, frontier_tails[(*iter)&1]);
    // these all reference the shared data structs
    size_t &curr_bin_index = shared_indexes[(*iter)&1];                          
    size_t &next_bin_index = shared_indexes[((*iter)+1)&1];                      
    size_t &curr_frontier_tail = frontier_tails[(*iter)&1];
    size_t &next_frontier_tail = frontier_tails[((*iter)+1)&1];
    //*aye = 0;
    //printf("PE %d | post init | iter = %d, curr_front = %lu\n", vp.pe, *iter, curr_frontier_tail);
    //printf("Pe %d | While start | iter = %d, shared_index[iter&1] = %lu\n", vp.pe, *iter, shared_indexes[(*iter)&1]);
    Partition<NodeID> fp(curr_frontier_tail);                                            // I don't know how to do a dynamic schedule nowait, so I'm just gonna naively partition the work
    //printf("Everyone here? PE %d\n", vp.pe);
    shmem_barrier_all();
    //printf("fp.start: %lu | fp.end: %lu\n", fp.start, fp.end);
    //printf("PE %d | Current frontier tail: %lu\n", vp.pe, frontier_tails[(*iter)&1]);
    for (size_t i = fp.start; i < fp.end; i++) {
      NodeID u = frontier[i];
      if (u >= vp.start && u < vp.end)
        udist = dist[u];
      else
        shmem_getmem(&udist, dist.begin()+vp.local_pos(u), sizeof(WeightT), vp.recv(u));                       // does this represent a synch issue? does it need to be made atomic somehow?
      //printf("Udist: %d\n", udist);
      //printf("dist[%d] = %d, curr_bin_index = %d, delta = %d\n", u, udist, static_cast<WeightT>(curr_bin_index), delta);
      if (udist >= delta * static_cast<WeightT>(curr_bin_index))
        RelaxEdges(g, u, delta, dist, local_bins, vp);
    }
    while (curr_bin_index < local_bins.size() &&
             !local_bins[curr_bin_index].empty() &&
             local_bins[curr_bin_index].size() < kBinSizeThreshold) {
      vector<NodeID> curr_bin_copy = local_bins[curr_bin_index];
      local_bins[curr_bin_index].resize(0);
      for (NodeID u : curr_bin_copy)
        RelaxEdges(g, u, delta, dist, local_bins, vp);
    }
    shmem_barrier_all();
    //shmem_global_exit(0);
    //exit(0);
    for (size_t i=curr_bin_index; i < local_bins.size(); i++) {         // this region could def be optimized, maybe each pe finds its local min b4 voting?
      if (!local_bins[i].empty()) {
        shmem_set_lock(VOTE_LOCK);                                            // critical region
        for (int p = 0; p < vp.npes; p++) {
          size_t temp = min(shared_indexes[((*iter)+1)&1], i);
      //    printf("PE %d is voting: temp = %lu\n", vp.pe, temp);
          shmem_putmem(shared_indexes+(((*iter)+1)&1), &temp, sizeof(size_t), p);
        }
        shmem_clear_lock(VOTE_LOCK);
        break;
      }      
    }
    shmem_barrier_all();
    //printf("PE %d | next chosen bin is %lu\n", vp.pe, shared_indexes[((*iter)+1)&1]);
    //shmem_global_exit(0);
    //exit(0);
    t.Stop();
    PrintStep(curr_bin_index, t.Millisecs(), curr_frontier_tail);
    t.Start();                                                          //REMEMBER! Everyone must particpate in timer and label printing stuff
    //if (shmem_test_lock(SINGLE_LOCK) == 0) {                                  // lock is not set, execute body 
    //curr_bin_index = kMaxBin;                                            // so since &curr_bin_index = shared_index, do changes to curr_bin_index affect (or should affect) the shared?
    shared_indexes[(*iter)&1] = kMaxBin;
    //curr_frontier_tail = 0;
    frontier_tails[(*iter)&1] = 0;
    //printf("PE %d has set curr bin to kmax and curr front to 0\n", vp.pe);
    //printf("PE %d | Curr bin index = %lu | curr_front = %lu | shared_indexes[iter&1] = %lu | shared_indexes[(iter+1)&1] = %lu\n", vp.pe, curr_bin_index, curr_frontier_tail, shared_indexes[(*iter)&1], shared_indexes[((*iter)+1)&1]);


    //shmem_clear_lock(SINGLE_LOCK);
    //}
    //printf("PE %d checking in\n", vp.pe);
    printf("PE %d | iter: %d \n", vp.pe, *iter);
    if (next_bin_index < local_bins.size()) {
      printf("Frontier contents b4 copy: ");
      for (WeightT t : frontier)
        printf("%d ", t);
      printf("\n");

      shmem_set_lock(SINGLE_LOCK);
      size_t copy_start;
      for (int i = 0; i < vp.npes; i++) {
        //if (vp.pe != i) {
        copy_start = shmem_ulong_atomic_fetch_add(&next_frontier_tail, local_bins[next_bin_index].size(), i);
        printf("Copy start: %lu\n", copy_start);
        printf("local.begin = %d\n", *(local_bins[next_bin_index].begin()));
        shmem_putmem(frontier.data() + copy_start, &*(local_bins[next_bin_index].begin()), sizeof(NodeID) * local_bins[next_bin_index].size(), i);
      }
      printf("Frontier contents after copy: ");
      for (WeightT t : frontier)
        printf("%d ", t);
      printf("\n");

      shmem_clear_lock(SINGLE_LOCK);
      local_bins[next_bin_index].resize(0);
    }

    printf("Pe %d | Prior to iter inc | iter = %d, shared_index[iter&1] = %lu, curr_front = %lu\n", vp.pe, *iter, shared_indexes[(*iter)&1], curr_frontier_tail);
    //iter++;
    if (vp.pe == 0) {
      for (int i = 0; i < vp.npes; i++) {
        shmem_int_atomic_inc(iter, i);
        printf("PE %d just incremented iter on PE %d\n", vp.pe, i);
      }
    }
    shmem_barrier_all();
    //printf("PE %d | While end | curr_front = %lu, next_front = %lu\n", vp.pe, 
    //printf("Pe %d | While end | iter = %d, shared_index[iter&1] = %lu, curr_front = %lu\n", vp.pe, *iter, shared_indexes[(*iter)&1], curr_frontier_tail);
    //printf("PE %d | While end | iter = %d, curr_frontier_tail = %lu\n", vp.pe, *iter, frontier_tails[(*iter)&1]); //curr_frontier_tail);
    printf("Pe %d | After iter inc | iter = %d, shared_index[iter&1] = %lu, curr_front = %lu\n", vp.pe, *iter, shared_indexes[(*iter)&1], curr_frontier_tail);
  }
  //printf("PE %d is here\n", vp.pe);
  if (shmem_test_lock(SINGLE_LOCK) == 0) {
    //std::cout << "took " << iter << " iterations" << endl;
    printf("PE %d | Took %d iterations\n", vp.pe, *iter);
    shmem_clear_lock(SINGLE_LOCK);
  }
  shmem_barrier_all();
  printf("end Pe %d\n", vp.pe);
  for (int i = vp.start; i < vp.end; i++)
    printf("PE %d | dist(%d) = %lu\n", vp.pe, i, dist[vp.local_pos(i)]);
  return dist;                                                          // only pe 0 has a correct distance array
}
 

/*pvector<WeightT> Shmem_DeltaStep(const WGraph &g, NodeID source, int delta, long* VOTE_LOCK, long* SINGLE_LOCK) {
  WeightT *old_dist, *new_dist;
  Timer t;
  pvector<WeightT> dist(g.num_nodes(), kDistInf, true);                                 // Initialize symmetric pvector with +inf - PE 0 maintains updated copy
  dist[source] = 0;
  pvector<NodeID> frontier(g.num_edges_directed(), true);                               // Symmetric frontier - only PE 0 updated (or maybe pe 0 has dist array, pe 1 has frontier)
  // two element arrays for double buffering curr=iter&1, next=(iter+1)&1
  size_t* shared_indexes = (size_t *) shmem_calloc(2, sizeof(size_t));                  // size_t shared_indexes[2] = {0, kMaxBin};
  shared_indexes[1] = kMaxBin;
  size_t* frontier_tails = (size_t *) shmem_calloc(2, sizeof(size_t));                   // size_t frontier_tails[2] = {1, 0};
  frontier_tails[0] = 1;
  frontier[0] = source;
  shmem_barrier_all();
  t.Start();
  vector<vector<NodeID> > local_bins(0);                                      // vector of vector of node ids (thread local?)
  size_t iter = 0;
  NodeID* i = (NodeID *) shmem_malloc(sizeof(NodeID));                                  // Shared counter for next node
  while (shared_indexes[iter&1] != kMaxBin) {                                 // if iter is odd, iter bitwise and (&) 1 is 1, if iter is even iter&1 is 0 and enter loop
    size_t &curr_bin_index = shared_indexes[iter&1];                          
    size_t &next_bin_index = shared_indexes[(iter+1)&1];                      
    size_t &curr_frontier_tail = frontier_tails[iter&1];
    size_t &next_frontier_tail = frontier_tails[(iter+1)&1];
    *i = 0;
    shmem_barrier_all();

    printf("curr_bin_index: %lu | next_bin_index: %lu | curr_front_tail: %lu | next_front_tail: %lu \n", curr_bin_index, next_bin_index, curr_frontier_tail, next_frontier_tail);
    while ((*i = shmem_int_atomic_fetch(i, 0) < curr_frontier_tail)) {
      int u = shmem_int_atomic_fetch(frontier.begin()+(*i), 0);                    // frontier only up to date on PE 0
      printf(" Pe %d thinks u = %d\n", shmem_my_pe(), u);
      shmem_global_exit(0);
      exit(0);
      shmem_int_get(dist.begin()+u, dist.begin()+u, 1, 0);                                      // PE 0 is the only PE with gauranteed up-to-date distance arrays 
      if (dist[u] >= delta * static_cast<int>(curr_bin_index)) {
        for (WNode wn : g.out_neigh(u)) {
          shmem_int_get(old_dist, dist.begin()+wn.v, 1, 0);
          shmem_int_get(new_dist, dist.begin()+u, 1, 0);
          *new_dist += wn.v;
          if (new_dist < old_dist) {
            bool changed_dist = true;
            while (!shmem_int_atomic_compare_swap(dist.begin()+wn.v, *old_dist, *new_dist, 0)) {       // While dist[wn.v] != old_dist, loop. when they become equal, dist[wn.v] is replaced with new_dist
              shmem_int_get(old_dist, dist.begin()+wn.v, 1, 0);                                             // use the local copy, or do we need PE 0's version?
              if (*old_dist <= *new_dist) {                                   
                changed_dist = false;
                break;
              }
            }
           if (changed_dist) {
                size_t dest_bin = (*new_dist)/delta;
                if (dest_bin >= local_bins.size()) {
                  local_bins.resize(dest_bin+1);
                }
                local_bins[dest_bin].push_back(wn.v);
              }
            }
          }
        }
        shmem_int_atomic_fetch_inc(i, 0);                                       // i++ on PE 0
      }
      shmem_barrier_all();                                                      // OpenMP synchs automatically after parallel for
      for (size_t i=curr_bin_index; i < local_bins.size(); i++) {
        if (!local_bins[i].empty()) {
          shmem_set_lock(VOTE_LOCK);                                            // critical region
          next_bin_index = min(next_bin_index, i);
          shmem_clear_lock(VOTE_LOCK);
          break;
        }
      }
      shmem_barrier_all();
      t.Stop();
      if (shmem_test_lock(SINGLE_LOCK) == 0) {                                  // lock is not set, execute body 
        PrintStep(curr_bin_index, t.Millisecs(), curr_frontier_tail);
        t.Start();
        curr_bin_index = kMaxBin;
        curr_frontier_tail = 0;
        shmem_clear_lock(SINGLE_LOCK);
      }
      if (next_bin_index < local_bins.size()) {
        //size_t copy_start = shmem_long_atomic_fetch_add(next_frontier_tail,
                      //                  local_bins[next_bin_index].size(), 0);
        //shmem_int_put(frontier.data() + copy_start, local_bins[next_bin_index].begin(), 
                    //  local_bins[next_bin_index].size(), 0);
        local_bins[next_bin_index].resize(0);
      }
      iter++;
      shmem_barrier_all();
    }
    if (shmem_test_lock(SINGLE_LOCK) == 0) {
      //o       std::cout << "took " << iter << " iterations" << endl;
      printf("Took %lu iterations\n", iter);
      shmem_clear_lock(SINGLE_LOCK);
    }
  return dist;                                                          // only pe 0 has a correct distance array
}*/

pvector<WeightT> DeltaStep(const WGraph &g, NodeID source, WeightT delta) {
  Timer t;
  pvector<WeightT> dist(g.num_nodes(), kDistInf);                          
  dist[source] = 0;
  pvector<NodeID> frontier(g.num_edges_directed());                               
  // two element arrays for double buffering curr=iter&1, next=(iter+1)&1
  size_t shared_indexes[2] = {0, kMaxBin};
  size_t frontier_tails[2] = {1, 0};
  frontier[0] = source;
  t.Start();
  #pragma omp parallel                                                          // all threads execute everything in this section
  {
    vector<vector<NodeID> > local_bins(0);                                      // vector of vector of node ids (thread local?)
    size_t iter = 0;
    while (shared_indexes[iter&1] != kMaxBin) {                                 // if iter is odd, iter bitwise and (&) 1 is 1, if iter is even iter&1 is 0
                                                                                // iter is even: shared[it&1] = 0, 0 != kMaxBin, condition is true
                                                                                // iter is odd: shared[it&1] = kMaxBin, condition is false
      size_t &curr_bin_index = shared_indexes[iter&1];                          
// If iter is even: &curr_bin_indx is a pointer to first elem of sharedindx, &next_bin_indx points to second elem, &curr_front points to first elem of frontier_tails, &next_front to second elem
      size_t &next_bin_index = shared_indexes[(iter+1)&1];                      
// If iter is odd: curr_bin_indx = kMaxBin, next_bin = 0, curr_front = 0, next_front = 1
      size_t &curr_frontier_tail = frontier_tails[iter&1];
      size_t &next_frontier_tail = frontier_tails[(iter+1)&1];
      #pragma omp for nowait schedule(dynamic, 64)                              // nowait: continue execution after exiting loop w/o synch
                                                                                // dynamic 64: assign iterations to threads in greedy manner, as soon as a thread is ready it takes the next iteration
      for (size_t i=0; i < curr_frontier_tail; i++) {
        NodeID u = frontier[i];
        if (dist[u] >= delta * static_cast<WeightT>(curr_bin_index)) {
          for (WNode wn : g.out_neigh(u)) {
            WeightT old_dist = dist[wn.v];
            WeightT new_dist = dist[u] + wn.w;
            if (new_dist < old_dist) {
              bool changed_dist = true;
              while (!compare_and_swap(dist[wn.v], old_dist, new_dist)) {       // While dist[wn.v] != old_dist, loop. when they become equal, dist[wn.v] is replaced with new_dist
                old_dist = dist[wn.v];                                          // maybe one pe should hold a master copy of the distances list that everyone pushes to and pulls from
                if (old_dist <= new_dist) {                                     // or divide the distances aray like the parent array?
                  changed_dist = false;
                  break;
                }
              }
              if (changed_dist) {
                size_t dest_bin = new_dist/delta;
                if (dest_bin >= local_bins.size()) {
                  local_bins.resize(dest_bin+1);
                }
                local_bins[dest_bin].push_back(wn.v);
              }
            }
          }
        }
      }
      for (size_t i=curr_bin_index; i < local_bins.size(); i++) {
        if (!local_bins[i].empty()) {
          #pragma omp critical                                          // changing next bin index changes the shared indexex (which is shared)
          next_bin_index = min(next_bin_index, i);
          break;
        }
      }
      #pragma omp barrier
      #pragma omp single nowait
      {
        t.Stop();
        PrintStep(curr_bin_index, t.Millisecs(), curr_frontier_tail);
        t.Start();
        curr_bin_index = kMaxBin;
        curr_frontier_tail = 0;
      }
      if (next_bin_index < local_bins.size()) {
        size_t copy_start = fetch_and_add(next_frontier_tail,
                                          local_bins[next_bin_index].size());
        copy(local_bins[next_bin_index].begin(),
             local_bins[next_bin_index].end(), frontier.data() + copy_start);
        local_bins[next_bin_index].resize(0);
      }
      iter++;
      #pragma omp barrier
    }
    #pragma omp single
    cout << "took " << iter << " iterations" << endl;
  }
  return dist;
}


void PrintSSSPStats(const WGraph &g, const pvector<WeightT> &dist) {
  auto NotInf = [](WeightT d) { return d != kDistInf; };
  int64_t num_reached = count_if(dist.begin(), dist.end(), NotInf);
  cout << "SSSP Tree reaches " << num_reached << " nodes" << endl;
}

// Compares against simple serial implementation
bool SSSPVerifier(const WGraph &g, NodeID source,
                  const pvector<WeightT> &dist_to_test) {
  // Serial Dijkstra implementation to get oracle distances
  pvector<WeightT> oracle_dist(g.num_nodes(), kDistInf);
  oracle_dist[source] = 0;
  typedef pair<WeightT, NodeID> WN;
  priority_queue<WN, vector<WN>, greater<WN>> mq;
  mq.push(make_pair(0, source));
  while (!mq.empty()) {
    WeightT td = mq.top().first;
    NodeID u = mq.top().second;
    mq.pop();
    if (td == oracle_dist[u]) {
      for (WNode wn : g.out_neigh(u)) {
        if (td + wn.w < oracle_dist[wn.v]) {
          oracle_dist[wn.v] = td + wn.w;
          mq.push(make_pair(td + wn.w, wn.v));
        }
      }
    }
  }
  // Report any mismatches
  bool all_ok = true;
  for (NodeID n : g.vertices()) {
    if (dist_to_test[n] != oracle_dist[n]) {
      cout << n << ": " << dist_to_test[n] << " != " << oracle_dist[n] << endl;
      all_ok = false;
    }
  }
  return all_ok;
}


int main(int argc, char* argv[]) {
  CLDelta<WeightT> cli(argc, argv, "single-source shortest-path");
  if (!cli.ParseArgs())
    return -1;

  char size_env[] = "SHMEM_SYMMETRIC_SIZE=1024M";
  putenv(size_env);

  static long PRINT_LOCK = 0;
  static long VOTE_LOCK = 0;
  static long SINGLE_LOCK = 0;

  shmem_init();

  static long pSync[SHMEM_REDUCE_SYNC_SIZE];
  static long pWrk[SHMEM_REDUCE_MIN_WRKDATA_SIZE];      
  static long long ll_pWrk[SHMEM_REDUCE_MIN_WRKDATA_SIZE];                 

  for (int i = 0; i < SHMEM_REDUCE_SYNC_SIZE; i++)
    pSync[i] = SHMEM_SYNC_VALUE;
  for (int i = 0; i < SHMEM_REDUCE_MIN_WRKDATA_SIZE; i++) {
    pWrk[i] = SHMEM_SYNC_VALUE;
    ll_pWrk[i] = SHMEM_SYNC_VALUE;
  }

  int pe = shmem_my_pe();
  int npes = shmem_n_pes();

  static long* PLOCKS = (long *) shmem_calloc(npes, sizeof(long));              // Access to shared resources controlled by a single pe is determined by a lock on each pe

  {
    WeightedBuilder b(cli);
    WGraph g = b.MakeGraph(pWrk, pSync);
    shmem_barrier_all();
    //g.PrintTopology(&PRINT_LOCK);
    SourcePicker<WGraph> sp(g, cli.start_vertex());
    /*Partition vp(g.num_nodes());
    printf("PE %d is init dist with width = %lu\n", vp.pe, vp.max_width);
    pvector<WeightT> dist(vp.max_width, kDistInf, true);                                 // Initialize symmetric partitioned pvector with +inf
    WeightT init = 0;
    NodeID source = 1;
    if (source >= vp.start && source < vp.end)
      dist[vp.local_pos(source)] = init;
    else
      shmem_putmem(dist.begin()+vp.local_pos(source), &init, sizeof(WeightT), vp.recv(source));
    printf("PE %d is init frontier with size = %lu\n", vp.pe, g.num_edges_directed());
    pvector<NodeID> frontier(g.num_edges_directed(), true);  
    if (shmem_int_atomic_compare_swap(dist.begin(), kDistInf, 15, 0))
      printf("swapped\n");*/
    auto SSSPBound = [&sp, &cli] (const WGraph &g) {
      return Shmem_DeltaStep(g, sp.PickNext(), cli.delta(), &VOTE_LOCK, &SINGLE_LOCK);
    };
    pvector<NodeID> parent = SSSPBound(g);
    shmem_barrier_all();
    printf("escape\n");  
  //SourcePicker<WGraph> vsp(g, cli.start_vertex());
  //auto VerifierBound = [&vsp] (const WGraph &g, const pvector<WeightT> &dist) {
    //return SSSPVerifier(g, vsp.PickNext(), dist);
  //};
  //BenchmarkKernel(cli, g, SSSPBound, PrintSSSPStats, VerifierBound);
  }
  shmem_finalize();
  return 0;
}





