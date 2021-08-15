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

/* 
DMM-GAPBS
Author: Zach Hansen
Adaptation Notes:
 - Work for top-down and bottom-up steps is divided between PES
 - All data structures except Bitmaps are distributed
*/

using namespace std;

// Assumes all PEs begin with synchronized front bitmaps, graph
// front is never updated within this function
// updates to parent arrays do not occur accross pe boundaries
// next bitmaps are synchronized at the end of the function
int64_t SHMEM_BUStep(const Graph &g, pvector<NodeID> &parent, Bitmap &front, Bitmap &next, int pe, int npes, long long *pWrk, long *psync, Partition<NodeID> vp) {
  int relative;
  next.reset();
  long long* awake_count = (long long *) shmem_malloc(sizeof(long long));               // Synchonization point
  *awake_count = 0;
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
int64_t SHMEM_TDStep(const Graph &g, pvector<NodeID> &parent, SlidingQueue<NodeID> &queue,
                      long *QLOCK, long *pSync, long long *pWrk, Partition<NodeID> vp) {
  long long* scout_count = (long long *) shmem_calloc(1, sizeof(long long));    // The global scout_count is in symmetric memory (calloc to init at 0)
  QueueBuffer<NodeID> lqueue(queue);                                            // Every PE maintains a queue buffer that updates the shared sliding queue
  for (NodeID u : queue) {                                                      // Each PE processes it's local portion of the active window
    NodeID curr_val;
    for (NodeID v : g.out_neigh(u)) {
      shmem_getmem(&curr_val, parent.begin()+vp.local_pos(v), sizeof(NodeID), vp.recv(v)); // v is the absolute location in the complete parent array
      if (curr_val < 0) {
        if (shmem_int_atomic_compare_swap(parent.begin()+vp.local_pos(v), curr_val, u, vp.recv(v)) == curr_val) {
          lqueue.push_back(v);
          *(scout_count) += -curr_val;
        }
      }
    }
  }
  lqueue.flush();
  shmem_longlong_sum_to_all(scout_count, scout_count, 1, 0, 0, vp.npes, pWrk, pSync);           // Reduction: + (represents a synchronization point)
  return (*scout_count);
}

// if the parent array is split accross many pes, it must be combined into one bitmap
// Bitmaps should be reset before this function
void QueueToBitmap(const SlidingQueue<NodeID> &queue, Bitmap &bm, long long* pWrk, long* pSync) {
  for (auto q_iter = queue.begin(); q_iter < queue.end(); q_iter++) {
    NodeID u = *q_iter;
    bm.set_bit(u);
  }
  shmem_barrier_all();
  bm.merge(pWrk, pSync);
  shmem_barrier_all();
}

// Assumes bitmaps are merged (synched) at function entry
void BitmapToQueue(const Graph &g, const Bitmap &bm, SlidingQueue<NodeID> &queue, long *QLOCK, int pe, int npes) {
  Partition<NodeID> vp(g.num_nodes());
  QueueBuffer<NodeID> lqueue(queue);
  for (NodeID n = vp.start; n < vp.end; n++) {
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
pvector<NodeID> InitParent(const Graph &g, NodeID source) {
  Partition<NodeID> p(g.num_nodes());
  pvector<NodeID> parent(p.max_width, true);                               // The last PE contains the remainding elements, so the symmetric parent array must be at least this large on each PE
  #pragma omp parallel for                                              // But even though the parent array is symmetric, the elements aren't the same across PEs
  for (NodeID n=p.start; n < p.end; n++)
    parent[p.local_pos(n)] = g.out_degree(n) != 0 ? -g.out_degree(n) : -1;
  if (source >= p.start && source < p.end)                                 // Source occurs in the local parent pvector
    parent[p.local_pos(source)] = source;
  return parent;
}

pvector<NodeID> DOBFS(const Graph &g, NodeID source, long *FRONTIER_LOCK, long long* pWrk, long* pSync, int alpha = 15, int beta = 18) {
  Partition<NodeID> vp(g.num_nodes());
  PrintStep("Source", static_cast<int64_t>(source));
  Timer t;
  t.Start();
  pvector<NodeID> parent = InitParent(g, source);
  t.Stop();
  PrintStep("i", t.Seconds());
  SlidingQueue<NodeID> frontier(vp.max_width, FRONTIER_LOCK);          // Partitioned symmetric queue
  if (shmem_my_pe() == 0)
    frontier.push_back(source);
  frontier.slide_window();
  Bitmap curr(g.num_nodes(), true);                     // Symmetric unpartitioned bitmap
  curr.reset();
  Bitmap front(g.num_nodes(), true);                    // Symmetric unpartitioned bitmap
  front.reset();                                        // All PEs are synched at this point
  int64_t edges_to_check = g.num_edges_directed();
  int64_t scout_count = g.out_degree(source);
  while (!frontier.empty()) {
    if (scout_count > edges_to_check / alpha) {
      int64_t awake_count, old_awake_count;
      TIME_OP(t, QueueToBitmap(frontier, front, pWrk, pSync));
      PrintStep("e", t.Seconds());
      awake_count = frontier.size();
      frontier.slide_window();
      do {
        t.Start();
        old_awake_count = awake_count;
        awake_count = SHMEM_BUStep(g, parent, front, curr, shmem_my_pe(), shmem_n_pes(), pWrk, pSync, vp);
        front.swap(curr);
        t.Stop();
        PrintStep("bu", t.Seconds(), awake_count);
      } while ((awake_count >= old_awake_count) || (awake_count > vp.max_width / beta));    // used to ac > g.num_nodes / beta, does vp.max_width make sense?
      TIME_OP(t, BitmapToQueue(g, front, frontier, FRONTIER_LOCK, shmem_my_pe(), shmem_n_pes()));
      PrintStep("c", t.Seconds());
      scout_count = 1;
    } else {
      t.Start();
      edges_to_check -= scout_count;
      scout_count = SHMEM_TDStep(g, parent, frontier, FRONTIER_LOCK, pSync, pWrk, vp);
      frontier.slide_window();
      t.Stop();
      PrintStep("td", t.Seconds(), frontier.size());
    }
  }
  return(parent);
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
  Partition<NodeID> vp(g.num_nodes());
  NodeID vertex;
  int* PRINTER = (int *) shmem_malloc(sizeof(int));
  *PRINTER = 0;
  shmem_barrier_all();
  shmem_int_wait_until(PRINTER, SHMEM_CMP_EQ, vp.pe);       // wait until previous PE puts your pe # in PRINTER
  ofstream shmem_out;
  shmem_out.open("bfs_output.txt", ios::app);
  for (NodeID n = vp.start; n < vp.end; n++) {
    vertex = parent[vp.local_pos(n)];
    if (vertex < -1) {
      shmem_out << -1 << endl;
    } else {
      shmem_out << parent[vp.local_pos(n)] << endl;
    }
  }
  shmem_out.close();
  if (!(vp.pe == vp.npes-1))
    shmem_int_p(PRINTER, vp.pe+1, vp.pe+1);
  shmem_free(PRINTER);
  return true;
}


int main(int argc, char* argv[]) {
  CLApp cli(argc, argv, "breadth-first search");
  if (!cli.ParseArgs())
    return -1;

  //char size_env[] = "SMA_SYMMETRIC_SIZE=16G";
  //putenv(size_env);

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
  // shmem_barrier_all();       for some ungodly reason i need a useless calloc before the next section or it deadlocks
  static long* PLOCKS = (long *) shmem_calloc(shmem_n_pes(), sizeof(long));

  {
    Builder b(cli, cli.do_verify());
    Graph g = b.MakeGraph(pWrk, pSync);
    shmem_barrier_all();
    SourcePicker<Graph> sp(g, cli.start_vertex());
    auto BFSBound = [&sp] (const Graph &g) { return DOBFS(g, sp.PickNext(), &FRONTIER_LOCK, ll_pWrk, pSync); };
    SourcePicker<Graph> vsp(g, cli.start_vertex());
    auto VerifierBound = [&vsp] (const Graph &g, const pvector<NodeID> &parent) {
      return BFSVerifier(g, vsp.PickNext(), parent);
    };
    BenchmarkKernel(cli, g, BFSBound, PrintBFSStats, VerifierBound);
  }                                                                                            // Extra scope to trigger deletion of graph, otherwise shmem destructor is screwy

  shmem_finalize();
  return 0;
}
