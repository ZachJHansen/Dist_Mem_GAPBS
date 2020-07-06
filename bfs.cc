// Copyright (c) 2015, The Regents of the University of California (Regents)
// See LICENSE.txt for license details

#include <iostream>
#include <vector>

#include "benchmark.h"
#include "bitmap.h"
#include "builder.h"
#include "command_line.h"
#include "graph.h"
#include "platform_atomics.h"
#include "pvector.h"
#include "sliding_queue.h"
#include "timer.h"


/*
GAP Benchmark Suite
Kernel: Breadth-First Search (BFS)
Author: Scott Beamer

Will return parent array for a BFS traversal from a source vertex

This BFS implementation makes use of the Direction-Optimizing approach [1].
It uses the alpha and beta parameters to determine whether to switch search
directions. For representing the frontier, it uses a SlidingQueue for the
top-down approach and a Bitmap for the bottom-up approach. To reduce
false-sharing for the top-down approach, thread-local QueueBuffer's are used.

To save time computing the number of edges exiting the frontier, this
implementation precomputes the degrees in bulk at the beginning by storing
them in parent array as negative numbers. Thus the encoding of parent is:
  parent[x] < 0 implies x is unvisited and parent[x] = -out_degree(x)
  parent[x] >= 0 implies x been visited

[1] Scott Beamer, Krste Asanović, and David Patterson. "Direction-Optimizing
    Breadth-First Search." International Conference on High Performance
    Computing, Networking, Storage and Analysis (SC), Salt Lake City, Utah,
    November 2012.
*/


using namespace std;

// Assumes all PEs begin with synchronized front bitmaps, graph
// front is never updated within this function
// updates to parent arrays do not occur accross pe boundaries
// next bitmaps are synchronized at the end of the function
int64_t SHMEM_BUStep(const Graph &g, pvector<NodeID> &parent, Bitmap &front, Bitmap &next, int pe, int npes, long long *pWrk, long *psync, Partition<NodeID> vp) {
  int relative;
  next.reset();
  long long* awake_count = (long long *) shmem_malloc(sizeof(long long));               // Synchonization point?
  *awake_count = 0;
  /*int parent_offset = g.num_nodes() / npes;
  int upper_bound, relative;                                                            // Each element of the complete parent array has a relative position in the local array
  int lower_bound = parent_offset * pe;                                                 // Distribute graph processing ~evenly
  if (pe == npes-1){
    upper_bound = g.num_nodes();                
  } else {
    upper_bound = lower_bound + parent_offset;
  }*/
  for (NodeID u = vp.start; u < vp.end; u++) {                                  // PE N has parent array[lower : upper] and is responsible for processing nodes lower-upper
    relative = vp.local_pos(u);
    if (parent[relative] < 0) {
      for (NodeID v : g.in_neigh(u)) {
        if (front.get_bit(v)) {
          parent[relative] = v;
          (*awake_count)++;
          next.set_bit(u);
          break;
        }
      }
    }
  }
  next.merge(pWrk, psync);                                                              // Synchronize local copies of bitmaps
  shmem_longlong_sum_to_all(awake_count, awake_count, 1, 0, 0, vp.npes, pWrk, psync);      // Reduction : +
  return(*awake_count);
}

// This assumes NodeIDs are integers, otherwise I'm not sure how to do the atomic compare and swap
// Assumes PLOCKS is an array of locks of length npes: so access to the parent array on each PE is controlled by a seperate lock
int64_t SHMEM_TDStep(const Graph &g, pvector<NodeID> &parent, SlidingQueue<NodeID> &queue, 
                      long *QLOCK, long *PLOCKS, long *pSync, long long *pWrk, Partition<NodeID> vp) {
  int pe = shmem_my_pe();
  int npes = shmem_n_pes();
  long long* scout_count = (long long *) shmem_calloc(1, sizeof(long long));    // The global scout_count is in symmetric memory (calloc to init at 0)
  QueueBuffer<NodeID> lqueue(queue);                                     // Every PE maintains a queue buffer that updates the shared sliding queue
  int queue_offset = queue.size() / npes;
  int parent_offset = g.num_nodes() / npes;
  int upper_bound, end, local_v, foreign_pe;
  int lower_bound = parent_offset * pe;                                         // Which members of the parent array are in the local parent array
  int start = queue_offset * pe;                                                // Divide processing of queue between PEs
  if (pe == npes-1){
    end = queue.size();
    upper_bound = g.num_nodes();                
  } else {
    end = start + queue_offset;
    upper_bound = lower_bound + parent_offset;
  }
  auto q_iter = queue.begin();
  q_iter += start;
  auto q_end = queue.begin();
  q_end += end;
  while (q_iter < q_end) {      // what if q size < npes?
    NodeID u = *q_iter;
    NodeID curr_val;
    for (NodeID v : g.out_neigh(u)) {
      if (v >= lower_bound && v < upper_bound) {                                        // The outgoing neighbor v of node u is in the local subset of the parent array
        shmem_set_lock(PLOCKS+pe);                                                      // A PE can lock itself to avoid simultaneous remote accesses to nodes in the local subset
        curr_val = parent[v-lower_bound];                                               // v is the absolute location in the complete parent array
        if (curr_val < 0) {
          parent[v-lower_bound] = u;                                                    // Update local subset with u
          lqueue.push_back(v);
          *(scout_count) += -curr_val;
        }  
        shmem_clear_lock(PLOCKS+pe);
      } else {                                                                          // v is in the parent array subset on a different PE
        foreign_pe = v / parent_offset;
        if (foreign_pe >= npes) {                                                       // The parent array on the last PE is > offset, thus foreign_pe could be >= npes
          foreign_pe = npes-1;
          local_v = v - parent_offset*foreign_pe;
        } else {
          local_v = v % parent_offset;
        }
        shmem_set_lock(PLOCKS+foreign_pe);                                             // Get exclusive access to the parent array on PE foreign_pe
        shmem_int_get(&curr_val, parent.begin()+local_v, 1, foreign_pe);                // Update curr_val with v from foreign pe's parent array
        if (curr_val < 0) {
          shmem_int_atomic_swap(parent.begin()+local_v, u, foreign_pe);                 // Update foregin pe with u
          lqueue.push_back(v);                                                          // The sliding queue and scouts get aggregated, so it shouldnt matter which PE updates them?
          *(scout_count) += -curr_val;
        }
        shmem_clear_lock(PLOCKS+foreign_pe);                                           
      }
    }
    q_iter++;
  }
  lqueue.flush();
  shmem_longlong_sum_to_all(scout_count, scout_count, 1, 0, 0, npes, pWrk, pSync);           // Reduction: + (represents a synchronization point)
  //printf("PE %d has scout count %lu\n", pe, *scout_count);
  return (*scout_count); 
}

// if the parent array is split accross many pes, it must be combined into one bitmap
// Bitmaps should be reset before this function
void QueueToBitmap(const SlidingQueue<NodeID> &queue, Bitmap &bm) {
  #pragma omp parallel for
  for (auto q_iter = queue.begin(); q_iter < queue.end(); q_iter++) {
    NodeID u = *q_iter;
    bm.set_bit_atomic(u);
  }
}

// Create afunction to calculate bounds, return them in a struct?
// Assumes bitmaps are merged (synched) at function entry
void BitmapToQueue(const Graph &g, const Bitmap &bm, SlidingQueue<NodeID> &queue, long *QLOCK, int pe, int npes) {
  int offset = g.num_nodes() / npes;
  int upper_bound;
  int lower_bound = offset * pe;                                                 // Distribute graph processing ~evenly
  if (pe == npes-1){
    upper_bound = g.num_nodes();                
  } else {
    upper_bound = lower_bound + offset;
  }
  QueueBuffer<NodeID> lqueue(queue);
  for (NodeID n = lower_bound; n < upper_bound; n++) {
      if (bm.get_bit(n)) {
        lqueue.push_back(n);
      }
  }
  lqueue.flush();
  queue.slide_window();                                                         // Slide window barrier_all PEs on function entry and exit 
}

// Partition parent array ~evenly across PEs (final PE gets remainder)
// Accessing node v on PE p means accessing node (n/k)*p + v in a complete parent array of n nodes and k PEs
// Similarly, node V in the complete parent array is the V%(n/k) element in the parent array of PE V/(n/k)
// (Unless pe = npes-1, then V is the V-(n/k)*p element in the npes-1 PE)
pvector<NodeID> InitParent(const Graph &g, NodeID source, int pe, int npes) {
  /*int start, end;
  int offset = g.num_nodes()/npes;
  size_t max_size = g.num_nodes() - (npes-1)*offset;
  start = offset * pe;
  if (pe == npes-1) {
    end = g.num_nodes();  
  } else {
    end = start + offset; 
  }*/
  Partition<NodeID> p(g.num_nodes());
  pvector<NodeID> parent(p.max_width, true);                               // The last PE contains the remainding elements, so the symmetric parent array must be at least this large on each PE
  #pragma omp parallel for                                              // But even though the parent array is symmetric, the elements aren't the same across PEs
  for (NodeID n=p.start; n < p.end; n++)
    parent[p.local_pos(n)] = g.out_degree(n) != 0 ? -g.out_degree(n) : -1;
  if (source >= p.start && source < p.end)                                 // Source occurs in the local parent pvector
    parent[p.local_pos(source)] = source;
  return parent;
}

pvector<NodeID> DOBFS(const Graph &g, NodeID source, long *FRONTIER_LOCK, long *PLOCKS, long long* pWrk, long* pSync, int alpha = 15, int beta = 18) {
  int pe = shmem_my_pe();
  int npes = shmem_n_pes();
  Partition<NodeID> vp(g.num_nodes());
  PrintStep("Source", static_cast<int64_t>(source));
  Timer t;
  t.Start();
  pvector<NodeID> parent = InitParent(g, source, pe, npes);
  t.Stop();
  PrintStep("i", t.Seconds());
  void* frontier_alloc = shmem_malloc(sizeof(SlidingQueue<NodeID>));                                
  SlidingQueue<NodeID>* frontier = new(frontier_alloc) SlidingQueue<NodeID>{(size_t) g.num_nodes(), FRONTIER_LOCK}; // too large?               
  if (pe == 0)
    frontier->push_back(source);
  frontier->slide_window();
  Bitmap curr(g.num_nodes(), true);                     // Symmetric bitmap
  curr.reset();
  Bitmap front(g.num_nodes(), true);                    // Symmetric bitmap
  front.reset();                                        // All PEs are synched at this point
  int64_t edges_to_check = g.num_edges_directed();
  int64_t scout_count = g.out_degree(source);
  //printf("Pe %d thin edges = %lu & scout = %lu\n", pe, edges_to_check, scout_count);
  /*shmem_barrier_all();
  shmem_global_exit(0);
  exit(0);*/
  while (!frontier->empty()) {
    if (scout_count > edges_to_check / alpha) {
      int64_t awake_count, old_awake_count;
      TIME_OP(t, QueueToBitmap(*frontier, front));
      PrintStep("e", t.Seconds());
      awake_count = frontier->size();
      frontier->slide_window();
      do {
        t.Start();
        old_awake_count = awake_count;
        awake_count = SHMEM_BUStep(g, parent, front, curr, pe, npes, pWrk, pSync, vp);
        front.swap(curr);
        t.Stop();
        PrintStep("bu", t.Seconds(), awake_count);
      } while ((awake_count >= old_awake_count) || (awake_count > vp.max_width / beta));    // used to ac > g.num_nodes / beta, does vp.max_width make sense?
      TIME_OP(t, BitmapToQueue(g, front, *frontier, FRONTIER_LOCK, pe, npes));
      PrintStep("c", t.Seconds());
      scout_count = 1;
    } else {
      t.Start();
      edges_to_check -= scout_count;
      scout_count = SHMEM_TDStep(g, parent, *frontier, FRONTIER_LOCK, PLOCKS, pSync, pWrk, vp);
      frontier->slide_window();
      t.Stop();
      PrintStep("td", t.Seconds(), frontier->size());
    }
  }
  //for (int i = 0; i < vp.max_width; i++) 
    //printf("PE %d | Node %d has parent %d\n", vp.pe, vp.global_pos(i), parent[i]);
  //for (p : parent)
    //printf("PE: %d | p - %d\n", vp.pe, p);
  // Do we care about the time required to combine parent arrays?
  // If not, it may be better to just fill a pvector with NodeIDs instead of 32 bit ints
  return(parent.combine(g.num_nodes(), pe, npes, pSync));
}

void PrintBFSStats(const Graph &g, const pvector<NodeID> &bfs_tree) {
  int64_t tree_size = 0;
  int64_t n_edges = 0;
  for (NodeID n : g.vertices()) {
    if (bfs_tree[n] >= 0) {
      n_edges += g.out_degree(n);
      tree_size++;
    }
  }
  cout << "BFS Tree has " << tree_size << " nodes and ";
  cout << n_edges << " edges" << endl;
}


// BFS verifier does a serial BFS from same source and asserts:
// - parent[source] = source
// - parent[v] = u  =>  depth[v] = depth[u] + 1 (except for source)
// - parent[v] = u  => there is edge from u to v
// - all vertices reachable from source have a parent

// BUT since we need a serial, non-partitioned verifier, we also need a complete
// graph built on PE 0. So verifier won't work on graphs that require partitioning
// output generated parent arrays to a file, then have the original BFS verify them
bool BFSVerifier(const Graph &g, NodeID source, const pvector<NodeID> &parent) {
  if (shmem_my_pe() == 0) {
    ofstream shmem_out;
    shmem_out.open("/home/zach/projects/Dist_Mem_GAPBS/Dist_Mem_GAPBS/shmem_output.txt", ios::app);
    for (auto it = parent.begin(); it < parent.end(); it++) {
      shmem_out << *it << endl;
  //  *it = 15;           // scramble
    }
  }
/*  pvector<int> depth(g.num_nodes(), -1);
  depth[source] = 0;
  vector<NodeID> to_visit;
  to_visit.reserve(g.num_nodes());
  to_visit.push_back(source);
  for (auto it = to_visit.begin(); it != to_visit.end(); it++) {
    NodeID u = *it;
    for (NodeID v : g.out_neigh(u)) {
      if (depth[v] == -1) {
        depth[v] = depth[u] + 1;
        to_visit.push_back(v);
      }
    }
  }
  for (NodeID u : g.vertices()) {
    if ((depth[u] != -1) && (parent[u] != -1)) {
      if (u == source) {
        if (!((parent[u] == u) && (depth[u] == 0))) {
          cout << "Source wrong" << endl;
          return false;
        }
        continue;
      }
      bool parent_found = false;
      for (NodeID v : g.in_neigh(u)) {
        if (v == parent[u]) {
          if (depth[v] != depth[u] - 1) {
            cout << "Wrong depths for " << u << " & " << v << endl;
            return false;
          }
          parent_found = true;
          break;
        }
      }
      if (!parent_found) {
        cout << "Couldn't find edge from " << parent[u] << " to " << u << endl;
        return false;
      }
    } else if (depth[u] != parent[u]) {
      cout << "Reachability mismatch" << endl;
      return false;
    }
  }*/
  return true;
}

int main(int argc, char* argv[]) {
  CLApp cli(argc, argv, "breadth-first search");
  if (!cli.ParseArgs())
    return -1;

  static long FRONTIER_LOCK = 0;                                                      // Create a mutex lock in symmetric memory to control access to the frontier
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

  char size_env[] = "SHMEM_SYMMETRIC_SIZE=40000000";
  putenv(size_env);
  int npes = shmem_n_pes();
  int pe = shmem_my_pe();
  static long* PLOCKS = (long *) shmem_calloc(shmem_n_pes(), sizeof(long));                                     // Access to shared resources controlled by a single pe is determined by a lock on each pe
  static long* LOCK = (long *) shmem_calloc(1, sizeof(long)); 
  {
    Builder b(cli);
    Graph g = b.MakeGraph(pWrk, pSync);
    //printf("Last check\n");
    //g.PrintTopology(LOCK);

    shmem_barrier_all();
    //printf("PE %d says out_index starts at %p\n", pe, g.out_index_);
    //printf("PE %d says out_neigh starts at %p\n", pe, g.out_neighbors_);
    /*int64_t x;
    for (v : g.vertices()) {
      x = g.out_degree(v);
      printf("PE %d | degree of %d: %lu\n", pe, v, x);
    }*/
    /*if (pe == 0) {
      for (int y = 0; y < 10; y++) {
        for (int x : g.in_neigh(y, 0))
          printf("%d neighbor: %d\n", y, x);
      }
    }*/
    SourcePicker<Graph> sp(g, cli.start_vertex());
    auto BFSBound = [&sp] (const Graph &g) { return DOBFS(g, sp.PickNext(), &FRONTIER_LOCK, PLOCKS, ll_pWrk, pSync); };
    SourcePicker<Graph> vsp(g, cli.start_vertex());
    auto VerifierBound = [&vsp] (const Graph &g, const pvector<NodeID> &parent) {
      return BFSVerifier(g, vsp.PickNext(), parent);
    }; 
    BenchmarkKernel(cli, g, BFSBound, PrintBFSStats, VerifierBound);
  }                                                                                            // Extra scope to trigger deletion of graph, otherwise shmem destructor is screwy
  shmem_free(PLOCKS);
  shmem_free(LOCK);
  shmem_finalize();
  return 0;
}


