// Copyright (c) 2015, The Regents of the University of California (Regents)
// See LICENSE.txt for license details

#include <functional>
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
#include "util.h"


/*
GAP Benchmark Suite
Kernel: Betweenness Centrality (BC)
Author: Scott Beamer

Will return array of approx betweenness centrality scores for each vertex

This BC implementation makes use of the Brandes [1] algorithm with
implementation optimizations from Madduri et al. [2]. It is only an approximate
because it does not compute the paths from every start vertex, but only a small
subset of them. Additionally, the scores are normalized to the range [0,1].

As an optimization to save memory, this implementation uses a Bitmap to hold
succ (list of successors) found during the BFS phase that are used in the back-
propagation phase.

[1] Ulrik Brandes. "A faster algorithm for betweenness centrality." Journal of
    Mathematical Sociology, 25(2):163–177, 2001.

[2] Kamesh Madduri, David Ediger, Karl Jiang, David A Bader, and Daniel
    Chavarria-Miranda. "A faster parallel algorithm and efficient multithreaded
    implementations for evaluating betweenness centrality on massive datasets."
    International Symposium on Parallel & Distributed Processing (IPDPS), 2009.
*/


// WARNING! MISSING ATOMIC CAS

using namespace std;
typedef float ScoreT;
typedef double CountT;


void PBFS(const Graph &g, NodeID source, pvector<CountT> &path_counts,
    Bitmap &succ, vector<SlidingQueue<NodeID>::iterator> &depth_index,
    SlidingQueue<NodeID> &queue) {
  pvector<NodeID> depths(g.num_nodes(), -1);
  depths[source] = 0;
  path_counts[source] = 1;
  queue.push_back(source);
  depth_index.push_back(queue.begin());
  queue.slide_window();
  const NodeID* g_out_start = g.out_neigh(0).begin();
  #pragma omp parallel
  {
    NodeID depth = 0;
    QueueBuffer<NodeID> lqueue(queue);
    while (!queue.empty()) {
      depth++;
      #pragma omp for schedule(dynamic, 64) nowait
      for (auto q_iter = queue.begin(); q_iter < queue.end(); q_iter++) {
        NodeID u = *q_iter;
        for (NodeID &v : g.out_neigh(u)) {
          if ((depths[v] == -1) &&
              (compare_and_swap(depths[v], static_cast<NodeID>(-1), depth))) {          // is the depths[v] == -1 comparison just for short circuit evaluation?
            lqueue.push_back(v);
          }
          if (depths[v] == depth) {
            succ.set_bit_atomic(&v - g_out_start);
            #pragma omp atomic
            path_counts[v] += path_counts[u];
          }
        }
      }
      lqueue.flush();
      #pragma omp barrier
      #pragma omp single
      {
        depth_index.push_back(queue.begin());
        queue.slide_window();
      }
    }
  }
  depth_index.push_back(queue.begin());
}

// path counts, depths are partitioned
void shmem_BFS(const Graph &g, NodeID source, pvector<CountT> &path_counts,
    Bitmap &succ, vector<SlidingQueue<NodeID>::iterator> &depth_index,
    SlidingQueue<NodeID> &queue, Partition<NodeID> vp, long* PATH_LOCK) {
  pvector<NodeID> depths(vp.max_width, -1);
  if (source >= vp.start && source < vp.end) {
    depths[vp.local_pos(source)] = 0;
    path_counts[vp.local_pos(source)] = 1;
    queue.push_back(source);                            // queue is symmetric but not partitioned, only one PE pushbacks
    depth_index.push_back(queue.begin());
  }
  queue.slide_window();                                 // synch point
  const NodeID* g_out_start = g.out_neigh(0).begin();   // this is the head of the neighbor array yeah? so how does partitioning affect this?
  NodeID depth = 0;
  QueueBuffer<NodeID> lqueue(queue);                    // each PE maintains a queue buffer 
  while (!queue.empty()) {
    depth++;                                            // Thread local so each PE should maintain and increment a copy
    Partition<NodeID> qp(queue.size());                 // Divide processing of queue - if npes > q.size then pe npes-1 handles entire queue, not optimal for large npes
    for (auto q_iter = qp.start; q_iter < qp.end; q_iter++) {
      NodeID u = *q_iter;
      for (NodeID &v : g.out_neigh(u)) {
        if (shmem_int_atomic_fetch(depths.begin()+vp.local_pos(v), vp.recv(v)) == -1) {                 // needs to be CAS 
          shmem_int_atomic_swap(depths.begin()+vp.local_pos(v), depth, vp.recv(v));
        //if ((depths[v] == -1) &&
        //    (compare_and_swap(depths[v], static_cast<NodeID>(-1), depth))) {          // is the depths[v] == -1 comparison just for short circuit evaluation?
          lqueue.push_back(v);
        }
        if (shmem_int_atomic_fetch(depths.begin()+vp.local_pos(v), vp.recv(v)) == depth) {
          //succ.set_bit_atomic(&v - g_out_start);                                // why does this have to be atomic? if its already set, resetting it shouldnt matter?
          succ.set_bit(&v - g_out_start);               // g_out might be screwed up by the partitioning...
          // this command is supposed to be atomic, but requires 2 shmem instructions at best:
          // shmem_atomic_fetch to get path_counts[u], and shmem_atomic_add to add that value to path_counts[v]
          // so is a lock the only solution?
          shmem_set_lock(PATH_LOCK);            // should be test lock since we dont want everyone to exec this sequentially
          CountT pc_u = shmem_double_g(path_counts.begin()+vp.local_pos(u), vp.recv(u));
          for (int i = 0; i < npes; i++)                // doesnt have to be a loop since only one PE has a copy of paths[v]
            shmem_double_atomic_add(path_counts.begin()+vp.local_pos(v), pc_u, i);
          shmem_clear_lock(PATH_LOCK);
          //path_counts[v] += path_counts[u];
        }
      }
    }
    lqueue.flush();
    if (vp.pe == 0)
      depth_index.push_back(queue.begin());
    queue.slide_window();                               // synch point
  }
  if (vp.pe == 0)
    depth_index.push_back(queue.begin());
}


pvector<ScoreT> Brandes(const Graph &g, SourcePicker<Graph> &sp,
                        NodeID num_iters, long* LOCK) {
  Timer t;
  t.Start();
  Partition<NodeID> vp(g.num_nodes());
  pvector<ScoreT> scores(vp.max_width, 0, true);                // symmetric partitioned pvector
  pvector<CountT> path_counts(vp.max_width, true);                   // symmetric partitioned pvector
  Bitmap succ(g.num_edges_directed());                          // symmetric non-partitioned bitmap
  vector<SlidingQueue<NodeID>::iterator> depth_index;
  SlidingQueue<NodeID> queue(g.num_nodes());                    // really wish i knew how to get away with a smaller queue. resize when necessary?
  t.Stop();
  PrintStep("a", t.Seconds());
  const NodeID* g_out_start = g.out_neigh(0).begin();
  for (NodeID iter=0; iter < num_iters; iter++) {
    NodeID source = sp.PickNext();
    cout << "source: " << source << endl;
    t.Start();
    path_counts.fill(0);
    depth_index.resize(0);
    queue.reset();
    succ.reset();
    PBFS(g, source, path_counts, succ, depth_index, queue, vp, LOCK);
    shmem_global_exit(0);
    exit(0);
    /*t.Stop();
    PrintStep("b", t.Seconds());
    pvector<ScoreT> deltas(g.num_nodes(), 0);
    t.Start();
    for (int d=depth_index.size()-2; d >= 0; d--) {
      #pragma omp parallel for schedule(dynamic, 64)
      for (auto it = depth_index[d]; it < depth_index[d+1]; it++) {
        NodeID u = *it;
        ScoreT delta_u = 0;
        for (NodeID &v : g.out_neigh(u)) {
          if (succ.get_bit(&v - g_out_start)) {
            delta_u += (path_counts[u] / path_counts[v]) * (1 + deltas[v]);
          }
        }
        deltas[u] = delta_u;
        scores[u] += delta_u;
      }
    }
    t.Stop();
    PrintStep("p", t.Seconds());*/
  }
  // normalize scores
  /*ScoreT biggest_score = 0;
  #pragma omp parallel for reduction(max : biggest_score)
  for (NodeID n=0; n < g.num_nodes(); n++)
    biggest_score = max(biggest_score, scores[n]);
  #pragma omp parallel for
  for (NodeID n=0; n < g.num_nodes(); n++)
    scores[n] = scores[n] / biggest_score;*/
  return scores;
}


void PrintTopScores(const Graph &g, const pvector<ScoreT> &scores) {
  vector<pair<NodeID, ScoreT>> score_pairs(g.num_nodes());
  for (NodeID n : g.vertices())
    score_pairs[n] = make_pair(n, scores[n]);
  int k = 5;
  vector<pair<ScoreT, NodeID>> top_k = TopK(score_pairs, k);
  for (auto kvp : top_k)
    cout << kvp.second << ":" << kvp.first << endl;
}


// Still uses Brandes algorithm, but has the following differences:
// - serial (no need for atomics or dynamic scheduling)
// - uses vector for BFS queue
// - regenerates farthest to closest traversal order from depths
// - regenerates successors from depths
bool BCVerifier(const Graph &g, SourcePicker<Graph> &sp, NodeID num_iters,
                const pvector<ScoreT> &scores_to_test) {
  pvector<ScoreT> scores(g.num_nodes(), 0);
  for (int iter=0; iter < num_iters; iter++) {
    NodeID source = sp.PickNext();
    // BFS phase, only records depth & path_counts
    pvector<int> depths(g.num_nodes(), -1);
    depths[source] = 0;
    vector<CountT> path_counts(g.num_nodes(), 0);
    path_counts[source] = 1;
    vector<NodeID> to_visit;
    to_visit.reserve(g.num_nodes());
    to_visit.push_back(source);
    for (auto it = to_visit.begin(); it != to_visit.end(); it++) {
      NodeID u = *it;
      for (NodeID v : g.out_neigh(u)) {
        if (depths[v] == -1) {
          depths[v] = depths[u] + 1;
          to_visit.push_back(v);
        }
        if (depths[v] == depths[u] + 1)
          path_counts[v] += path_counts[u];
      }
    }
    // Get lists of vertices at each depth
    vector<vector<NodeID>> verts_at_depth;
    for (NodeID n : g.vertices()) {
      if (depths[n] != -1) {
        if (depths[n] >= static_cast<int>(verts_at_depth.size()))
          verts_at_depth.resize(depths[n] + 1);
        verts_at_depth[depths[n]].push_back(n);
      }
    }
    // Going from farthest to clostest, compute "depencies" (deltas)
    pvector<ScoreT> deltas(g.num_nodes(), 0);
    for (int depth=verts_at_depth.size()-1; depth >= 0; depth--) {
      for (NodeID u : verts_at_depth[depth]) {
        for (NodeID v : g.out_neigh(u)) {
          if (depths[v] == depths[u] + 1) {
            deltas[u] += (path_counts[u] / path_counts[v]) * (1 + deltas[v]);
          }
        }
        scores[u] += deltas[u];
      }
    }
  }
  // Normalize scores
  ScoreT biggest_score = *max_element(scores.begin(), scores.end());
  for (NodeID n : g.vertices())
    scores[n] = scores[n] / biggest_score;
  // Compare scores
  bool all_ok = true;
  for (NodeID n : g.vertices()) {
    ScoreT delta = abs(scores_to_test[n] - scores[n]);
    if (delta > std::numeric_limits<ScoreT>::epsilon()) {
      cout << n << ": " << scores[n] << " != " << scores_to_test[n];
      cout << "(" << delta << ")" << endl;
      all_ok = false;
    }
  }
  return all_ok;
}


int main(int argc, char* argv[]) {
  CLIterApp cli(argc, argv, "betweenness-centrality", 1);
  if (!cli.ParseArgs())
    return -1;

  if (cli.num_iters() > 1 && cli.start_vertex() != -1)
    cout << "Warning: iterating from same source (-r & -i)" << endl;

  shmem_init();

  long* LOCK = (long*) shmem_calloc(1, sizeof(long));
  static long pSync[SHMEM_REDUCE_SYNC_SIZE];
  static long lng_pWrk[SHMEM_REDUCE_MIN_WRKDATA_SIZE];      
  static double dbl_pWrk[SHMEM_REDUCE_MIN_WRKDATA_SIZE];                 

  for (int i = 0; i < SHMEM_REDUCE_SYNC_SIZE; i++)
    pSync[i] = SHMEM_SYNC_VALUE;
  for (int i = 0; i < SHMEM_REDUCE_MIN_WRKDATA_SIZE; i++) {
    lng_pWrk[i] = SHMEM_SYNC_VALUE;
    dbl_pWrk[i] = SHMEM_SYNC_VALUE;
  }

  {
    Builder b(cli);
    Graph g = b.MakeGraph(pSync, lng_pWrk);
    SourcePicker<Graph> sp(g, cli.start_vertex());
    auto BCBound =
      [&sp, &cli] (const Graph &g) { return Brandes(g, sp, cli.num_iters(), LOCK); };
    SourcePicker<Graph> vsp(g, cli.start_vertex());
    auto VerifierBound = [&vsp, &cli] (const Graph &g,
                                     const pvector<ScoreT> &scores) {
      return BCVerifier(g, vsp, cli.num_iters(), scores);
    };
    BenchmarkKernel(cli, g, BCBound, PrintTopScores, VerifierBound);
    shmem_free(LOCK);
  }
  shmem_finalize();
  return 0;
}
