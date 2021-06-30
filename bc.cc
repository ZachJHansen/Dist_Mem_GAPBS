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


using namespace std;
typedef float ScoreT;
typedef unsigned long long CountT;


// path counts, depths, queue are partitioned
void shmem_BFS(const Graph &g, NodeID source, pvector<CountT> &path_counts,
    Bitmap &succ, vector<SlidingQueue<NodeID>::iterator> &depth_index,
    SlidingQueue<NodeID> &queue, Partition<NodeID> vp) {

  static long BFSpSync[SHMEM_REDUCE_SYNC_SIZE];
  static long long BFSpWrk[SHMEM_REDUCE_MIN_WRKDATA_SIZE];
  for (int i = 0; i < SHMEM_REDUCE_SYNC_SIZE; i++)
    BFSpSync[i] = SHMEM_SYNC_VALUE;
  for (int i = 0; i < SHMEM_REDUCE_MIN_WRKDATA_SIZE; i++)
    BFSpWrk[i] = SHMEM_SYNC_VALUE;

  pvector<long> depths(vp.max_width, -1, true);
  if (source >= vp.start && source < vp.end) {
    depths[vp.local_pos(source)] = 0;
    path_counts[vp.local_pos(source)] = 1;
    queue.push_back(source);                            // queue is symmetric and partitioned
  }
  depth_index.push_back(queue.begin());
  queue.slide_window();                                 // synch point
  long depth = 0;
  CountT old_v;
  QueueBuffer<NodeID> lqueue(queue);                    // each PE maintains a queue buffer 
  shmem_barrier_all();
  while (!queue.empty()) {                              // while at least one PE has a non-empty frontier
    depth++;                                            // Thread local so each PE should maintain and increment a copy
    for (NodeID u : queue) {                            // Each PE processes their entire active portion of the queue
      NodeID neighbor_counter = 0;
      for (NodeID v : g.out_neigh(u)) {
        NodeID lp_v = vp.local_pos(v);
        if ((shmem_long_atomic_compare_swap(depths.begin()+lp_v, static_cast<long>(-1), depth, vp.recv(v)) == -1)) {
          lqueue.push_back(v);
        }
        if (shmem_long_g(depths.begin()+lp_v, vp.recv(v)) == depth) {
          succ.set_bit_partitioned(g, u, neighbor_counter, vp);
          old_v = shmem_ulonglong_atomic_fetch_add(path_counts.begin()+lp_v, shmem_ulonglong_atomic_fetch(path_counts.begin()+vp.local_pos(u), vp.recv(u)), vp.recv(v));
        }
        neighbor_counter++;
      }
    }
    lqueue.flush();
    shmem_barrier_all();
    depth_index.push_back(queue.begin());
    queue.slide_window();                               // synch point
  }
  shmem_barrier_all();
  succ.merge(BFSpWrk, BFSpSync);
  depth_index.push_back(queue.begin());
}


pvector<ScoreT> Brandes(const Graph &g, SourcePicker<Graph> &sp,
                        NodeID num_iters, long* QLOCK, float* pWrk, long* pSync) {
  Timer t;
  t.Start();
  Partition<NodeID> vp(g.num_nodes());
  pvector<ScoreT> scores(vp.max_width, 0, true);                        // symmetric partitioned pvector
  pvector<CountT> path_counts(vp.max_width, true);                      // symmetric partitioned pvector
  Bitmap succ(g.num_edges_directed(), true);                            // symmetric non-partitioned bitmap
  succ.init_bitmap_offsets(g, vp);
  vector<SlidingQueue<NodeID>::iterator> depth_index;
  SlidingQueue<NodeID> queue(vp.max_width, QLOCK);                      // symmetric partitioned frontier (queue) 
  t.Stop();
  PrintStep("a", t.Seconds());
  for (NodeID iter=0; iter < num_iters; iter++) {
    NodeID source = sp.PickNext();
    cout << "source: " << source << endl;
    t.Start();
    path_counts.fill(0);
    depth_index.resize(0);
    queue.reset();
    succ.reset();
    shmem_BFS(g, source, path_counts, succ, depth_index, queue, vp);
    t.Stop();
    PrintStep("b", t.Seconds());
    pvector<ScoreT> deltas(vp.max_width, 0, true);                      // symmetric partitioned pvector 
    t.Start();
    for (long d=depth_index.size()-2; d >= 0; d--) {                    // the partitioned frontier is divided into regions based on depth
      ScoreT delta_u;
      NodeID u, lp_u, lp_v;
      for (auto it = depth_index[d]; it < depth_index[d+1]; it++) {     // process all nodes in local copy of frontier within given depth
        u = *it;
        lp_u = vp.local_pos(u);
        delta_u = 0;
        NodeID neighbor_counter = 0;
        for (NodeID v : g.out_neigh(u)) {
          lp_v = vp.local_pos(v);
          if (succ.get_bit_partitioned(g, u, neighbor_counter, vp)) {
            double pc_v = (double) shmem_ulonglong_g(path_counts.begin()+lp_v, vp.recv(v));
            delta_u += (((double) shmem_ulonglong_g(path_counts.begin()+lp_u, vp.recv(u))) / pc_v) * (1 + shmem_float_g(deltas.begin()+lp_v, vp.recv(v)));
            }
            neighbor_counter++;
          }
          if (vp.pe == vp.recv(u)) {
            deltas[lp_u] = delta_u;
            scores[lp_u] += delta_u;
          } else {
            shmem_float_p(deltas.begin()+lp_u, delta_u, vp.recv(u));
            shmem_float_p(scores.begin()+lp_u, shmem_float_g(scores.begin()+lp_u, vp.recv(u))+delta_u, vp.recv(u));
          }
        }
        shmem_barrier_all();                                              // synchronize between depths
      }
    t.Stop();
    PrintStep("p", t.Seconds());
  }

  ScoreT* biggest_score = (float *) shmem_malloc(sizeof(float));        // Normalize scores
  *biggest_score = 0;
  for (NodeID n=vp.start; n < vp.end; n++)
    *biggest_score = max(*biggest_score, scores[vp.local_pos(n)]);
  shmem_barrier_all();
  shmem_float_max_to_all(biggest_score, biggest_score, 1, 0, 0, vp.npes, pWrk, pSync);
  for (NodeID n=vp.start; n < vp.end; n++)
    scores[vp.local_pos(n)] = scores[vp.local_pos(n)] / *biggest_score;
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


// Prints result to file to be read by original verifier
bool BCVerifier(const Graph &g, SourcePicker<Graph> &sp, NodeID num_iters,
                const pvector<ScoreT> &scores_to_test) {
  Partition<NodeID> vp(g.num_nodes());
  int* PRINTER = (int *) shmem_malloc(sizeof(int));
  *PRINTER = 0;
  shmem_barrier_all();
  shmem_int_wait_until(PRINTER, SHMEM_CMP_EQ, vp.pe);       // wait until previous PE puts your pe # in PRINTER
  ofstream shmem_out;
  shmem_out.precision(17);
  shmem_out.open("bc_output.txt", ios::app);
  for (NodeID n = vp.start; n < vp.end; n++) {
    shmem_out << scores_to_test[vp.local_pos(n)] << endl;
  }
  shmem_out.close();
  if (!(vp.pe == vp.npes-1))
    shmem_int_p(PRINTER, vp.pe+1, vp.pe+1);
  shmem_free(PRINTER);
  return true;
}


int main(int argc, char* argv[]) {
  CLIterApp cli(argc, argv, "betweenness-centrality", 1);
  if (!cli.ParseArgs())
    return -1;

  if (cli.num_iters() > 1 && cli.start_vertex() != -1)
    cout << "Warning: iterating from same source (-r & -i)" << endl;

  char size_env[] = "SMA_SYMMETRIC_SIZE=16G";
  putenv(size_env);

  static long QLOCK = 0;

  shmem_init();

  static long pSync[SHMEM_REDUCE_SYNC_SIZE];
  static long lng_pWrk[SHMEM_REDUCE_MIN_WRKDATA_SIZE];
  static float flt_pWrk[SHMEM_REDUCE_MIN_WRKDATA_SIZE];

  for (int i = 0; i < SHMEM_REDUCE_SYNC_SIZE; i++)
    pSync[i] = SHMEM_SYNC_VALUE;
  for (int i = 0; i < SHMEM_REDUCE_MIN_WRKDATA_SIZE; i++) {
    lng_pWrk[i] = SHMEM_SYNC_VALUE;
    flt_pWrk[i] = SHMEM_SYNC_VALUE;
  }

  {
    Builder b(cli, cli.do_verify());
    Graph g = b.MakeGraph(lng_pWrk, pSync);
    shmem_barrier_all();
    SourcePicker<Graph> sp(g, cli.start_vertex());
    auto BCBound =
      [&sp, &cli] (const Graph &g) { return Brandes(g, sp, cli.num_iters(), &QLOCK, flt_pWrk, pSync); };
    SourcePicker<Graph> vsp(g, cli.start_vertex());
    auto VerifierBound = [&vsp, &cli] (const Graph &g,
                                     const pvector<ScoreT> &scores) {
      return BCVerifier(g, vsp, cli.num_iters(), scores);
    };
    BenchmarkKernel(cli, g, BCBound, PrintTopScores, VerifierBound);
  }
  shmem_finalize();
  return 0;
}
