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

int64_t BUStep(const Graph &g, pvector<NodeID> &parent, Bitmap &front,
               Bitmap &next) {
  int64_t awake_count = 0;
  next.reset();
  #pragma omp parallel for reduction(+ : awake_count) schedule(dynamic, 1024)
  for (NodeID u=0; u < g.num_nodes(); u++) {
    if (parent[u] < 0) {
      for (NodeID v : g.in_neigh(u)) {
        if (front.get_bit(v)) {
          parent[u] = v;
          awake_count++;
          next.set_bit(u);
          break;
        }
      }
    }
  }
  return awake_count;
}

// This assumes NodeIDs are integers, otherwise I'm not sure how to do the atomic compare and swap
// Assumes PLOCKS is an array of locks of length npes: so access to the parent array on each PE is controlled by a seperate lock
int64_t SHMEM_TDStep(const Graph &g, pvector<NodeID> &parent, SlidingQueue<NodeID> &queue, 
                      long *QLOCK, long *PLOCKS, long *pSync, int *pwrk) {
  int pe = shmem_my_pe();
  int npes = shmem_n_pes();
  int64_t local_scout;                                                          // Every PE maintains a scout count for their subset of the parent array 
  int* scout_count = (int *) shmem_malloc(sizeof(int));                         // The global scout_count is in symmetric memory
  QueueBuffer<NodeID> lqueue(queue, QLOCK);                                     // Every PE maintains a queue buffer that updates the shared sliding queue
  int queue_offset = queue.size() / npes;
  int parent_offset = g.num_nodes() / npes;
  int upper_bound, end;
  int lower_bound = parent_offset * pe;                                         // Which members of the parent array are in the local parent array
  int start = queue_offset * pe;                                                // Divide processing of queue between PEs
  if (pe == npes-1){
    end = queue.size();
    upper_bound = g.num_nodes();                
  } else {
    end = start + offset;
    upper_bound = lower_bound + parent_offset;
  }
  auto q_iter = queue.begin();
  q_iter += start;
  auto q_end = queue.begin();
  q_end += end;
  printf("PE %d | Start: %d | End: %d\n", pe, start, end);
  while (q_iter < q_end) {
    NodeID u = *q_iter;
    //printf("U: %d\n", u);
    for (NodeID v : g.out_neigh(u)) {
      if (v >= lower_bound && v < upper_bound) {                                        // The outgoing neighbor v of node u is in the local subset of the parent array
        NodeID curr_val = parent[v-lower_bound];
        if (curr_val < 0) {
          if (shmem_int_atomic_compare_swap(parent[v-lower_bound], curr_val, u, pe)) {  // If a remote put has not updated parent[v], then replace value of parent[v] with u
            lqueue.push_back(v);
            local_scout += -curr_val;
          }  
        }
      } else {                                                                          // v is in the parent array subset on a different PE
        NodeID curr_val;
        int local_v;
        int foreign_pe = v / parent_offset;
        if (foreign_pe >= npes) {                                                       // The parent array on the last PE is > offset, thus foreign_pe could be >= npes
          foreign_pe = npes-1;
          local_v = v - parent_offset*foreign_pe;
        } else {
          local_v = v % parent_offset;
        }
        shmem_set_lock(PLOCKS[foreign_pe]);                                             // Get exclusive access to the parent array on PE foreign_pe
        shmem_int_get(curr_val, parent+local_v, 1, foreign_pe); 
        if (curr_val < 0) {
          shmem_int_atomic_swap(parent[local_v], u, foreign_pe);
          lqueue.push_back(v);                                                          // The sliding queue and scouts get aggregated, so it shouldnt matter which PE updates them?
          local_scout += -curr_val;
        }
        shmem_clear_lock(PLOCKS[foreign_pe]);                                           
      }
    }
    q_iter++;
  }
  lqueue.flush();
  shmem_int_sum_to_all(scout_count, local_scout, 1, 0, 0, npes, pwrk, pSync);           // Reduction: + (represents a synchronization point)
  printf("PE: %d | Local: %d | Total: %d\n", pe, *local_scout_count, *scout_count);
  return ((int64_t) *scout_count); 
}
 
int64_t TDStep(const Graph &g, pvector<NodeID> &parent, SlidingQueue<NodeID> &queue, long *QLOCK) {
  int64_t scout_count = 0;
  #pragma omp parallel
  {
    QueueBuffer<NodeID> lqueue(queue, QLOCK);
    #pragma omp for reduction(+ : scout_count)
    for (auto q_iter = queue.begin(); q_iter < queue.end(); q_iter++) {
      NodeID u = *q_iter;
      for (NodeID v : g.out_neigh(u)) {
        NodeID curr_val = parent[v];
        if (curr_val < 0) {
          if (compare_and_swap(parent[v], curr_val, u)) {
            lqueue.push_back(v);
            scout_count += -curr_val;
          }
        }
      }
    }
    lqueue.flush();
  }
  return scout_count;
}


void QueueToBitmap(const SlidingQueue<NodeID> &queue, Bitmap &bm) {
  #pragma omp parallel for
  for (auto q_iter = queue.begin(); q_iter < queue.end(); q_iter++) {
    NodeID u = *q_iter;
    bm.set_bit_atomic(u);
  }
}

void BitmapToQueue(const Graph &g, const Bitmap &bm, SlidingQueue<NodeID> &queue, long *QLOCK) {
  #pragma omp parallel
  {
    QueueBuffer<NodeID> lqueue(queue, QLOCK);
    #pragma omp for
    for (NodeID n=0; n < g.num_nodes(); n++)
      if (bm.get_bit(n))
        lqueue.push_back(n);
    lqueue.flush();
  }
  queue.slide_window();
}

// Partition parent array ~evenly across PEs (final PE gets remainder)
// Accessing node v on PE p means accessing node (n/k)*p + v in a complete parent array of n nodes and k PEs
// Similarly, node V in the complete parent array is the V%(n/k) element in the parent array of PE V/(n/k)
// (Unless pe = npes-1, then V is the V-(n/k)*p element in the npes-1 PE)
pvector<NodeID> InitParent(const Graph &g) {
  int start, end;
  int pe = shmem_my_pe();
  int npes = shmem_n_pes();
  int offset = g.num_nodes()/npes;
  size_t max_size = g.num_nodes() - (npes-1)*offset;
  start = offset * pe;
  if (pe == npes-1) {
    end = g.num_nodes();  
  } else {
    end = start + offset; 
  }
  pvector<NodeID> parent(max_size, true);                               // The last PE contains the remainding elements, so the symmetric parent array must be at least this large on each PE
  #pragma omp parallel for                                              // But even though the parent array is symmetric, the elements aren't the same across PEs
  for (NodeID n=start; n < end; n++)
    parent[n-start] = g.out_degree(n) != 0 ? -g.out_degree(n) : -1;
  return parent;
}

pvector<NodeID> DOBFS(const Graph &g, NodeID source, long *FRONTIER_LOCK, 
                      long *pSync, int *pwrk, int alpha = 15, int beta = 18) {
  PrintStep("Source", static_cast<int64_t>(source));
  Timer t;
  t.Start();
  pvector<NodeID> parent = InitParent(g);
  static long* PLOCKS = (long *) shmem_calloc(shmem_n_pes(), sizeof(long));
  t.Stop();
  PrintStep("i", t.Seconds());
  parent[source] = source;
  void* frontier_alloc = shmem_malloc(sizeof(SlidingQueue<NodeID>));                                
  SlidingQueue<NodeID>* frontier = new(frontier_alloc) SlidingQueue<NodeID>{(size_t) g.num_nodes()};                                  
//  SlidingQueue<NodeID> queue(g.num_nodes());
  frontier->push_back(source);
  frontier->slide_window();
  Bitmap curr(g.num_nodes());
  curr.reset();
  Bitmap front(g.num_nodes());
  front.reset();
  int64_t edges_to_check = g.num_edges_directed();
  int64_t scout_count = g.out_degree(source);
  printf("Initial edges to check: %d | Initial scout: %d\n", edges_to_check, scout_count);
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
        awake_count = BUStep(g, parent, front, curr);
        front.swap(curr);
        t.Stop();
        PrintStep("bu", t.Seconds(), awake_count);
      } while ((awake_count >= old_awake_count) ||
               (awake_count > g.num_nodes() / beta));
      TIME_OP(t, BitmapToQueue(g, front, *frontier, FRONTIER_LOCK));
      PrintStep("c", t.Seconds());
      scout_count = 1;
    } else {
      t.Start();
      edges_to_check -= scout_count;
      printf("Beginning TDStep\n");
      //scout_count = TDStep(g, parent, *frontier, FRONTIER_LOCK);   
      scout_count = SHMEM_TDStep(g, parent, *frontier, FRONTIER_LOCK, PLOCKS, pSync, pwrk);
      //exit(0);
      frontier->slide_window();
      t.Stop();
      PrintStep("td", t.Seconds(), frontier->size());
    }
  }
  #pragma omp parallel for
  for (NodeID n = 0; n < g.num_nodes(); n++)
    if (parent[n] < -1)
      parent[n] = -1;
  return parent;
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
bool BFSVerifier(const Graph &g, NodeID source,
                 const pvector<NodeID> &parent) {
  pvector<int> depth(g.num_nodes(), -1);
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
  }
  return true;
}


int main(int argc, char* argv[]) {
  CLApp cli(argc, argv, "breadth-first search");
  if (!cli.ParseArgs())
    return -1;

  static long FRONTIER_LOCK = 0;                                                      // Create a mutex lock in symmetric memory to control access to the frontier
  static int pwrk[SHMEM_REDUCE_MIN_WRKDATA_SIZE];
  static long pSync[SHMEM_REDUCE_SYNC_SIZE];
  for (int i = 0; i < SHMEM_REDUCE_SYNC_SIZE; i++)
    pSync[i] = SHMEM_SYNC_VALUE;
  for (int i = 0; i < SHMEM_REDUCE_MIN_WRKDATA_SIZE; i++)
    pwrk[i] = SHMEM_SYNC_VALUE;
  /*static long pSync[SHMEM_COLLECT_SYNC_SIZE];
  for (int i = 0; i < SHMEM_COLLECT_SYNC_SIZE; i++)
    pSync[i] = SHMEM_SYNC_VALUE;    
  */
  shmem_init();

  int npes = shmem_n_pes();
  int pe = shmem_my_pe();

  {
    void* builder_alloc = shmem_malloc(sizeof(Builder));
    Builder* b = new(builder_alloc) Builder{cli};
    Graph g = b->MakeGraph();
    g.PrintTopology();
    shmem_barrier_all();
    printf("Directed: %d | Num Edges: %d | Num Directed Edges: %d\n", g.directed(), g.num_edges(), g.num_edges_directed());    
    SourcePicker<Graph> sp(g, cli.start_vertex());
  
   /* void* frontier_alloc = shmem_malloc(sizeof(SlidingQueue<NodeID>));                                
    SlidingQueue<NodeID>* frontier = new(frontier_alloc) SlidingQueue<NodeID>{(size_t) g.num_nodes()};                                  
    for (int i = 0; i < g.num_nodes(); i++)
      frontier->push_back(i+10);
    frontier->slide_window();
    int64_t thing = SHMEM_TDStep(g, parent, *frontier, &FRONTIER_LOCK, PLOCKS);
    printf("Thing: %d\n", thing);*/
    auto BFSBound = [&sp] (const Graph &g) { return DOBFS(g, sp.PickNext(), &FRONTIER_LOCK, pSync, pwrk); };
    BFSBound(g);
   /*  SourcePicker<Graph> vsp(g, cli.start_vertex());
    auto VerifierBound = [&vsp] (const Graph &g, const pvector<NodeID> &parent) {
    return BFSVerifier(g, vsp.PickNext(), parent);
    };
    BenchmarkKernel(cli, g, BFSBound, PrintBFSStats, VerifierBound);*/
    shmem_free(b);
  }                                                                                            // Extra scope to trigger deletion of graph, otherwise shmem destructor is screwy
  shmem_finalize();
  return 0;
}


