// Copyright (c) 2015, The Regents of the University of California (Regents)
// See LICENSE.txt for license details

#include <algorithm>
#include <iostream>
#include <vector>
#include <iomanip>

#include "benchmark.h"
#include "builder.h"
#include "command_line.h"
#include "graph.h"
#include "pvector.h"

/*
GAP Benchmark Suite
Kernel: PageRank (PR)
Author: Scott Beamer

Will return pagerank scores for all vertices once total change < epsilon

This PR implementation uses the traditional iterative approach. This is done
to ease comparisons to other implementations (often use same algorithm), but
it is not necesarily the fastest way to implement it. It does perform the
updates in the pull direction to remove the need for atomics.
*/

/*
DMM-GAPBS
Author: Zach Hansen
Adaptation Notes:
 - Naively divide the work for each iteration
 - Precision issues are still making verification difficult
*/

using namespace std;

typedef float ScoreT;
const float kDamp = 0.85;

pvector<ScoreT> PageRankPull(const Graph &g, int max_iters,
                             long* pSync, double* pWrk, double epsilon = 0) {
  ScoreT temp, incoming_total, old_score;
  double* error = (double *) shmem_malloc(sizeof(double));
  const ScoreT init_score = 1.0f / g.num_nodes();
  const ScoreT base_score = (1.0f - kDamp) / g.num_nodes();
  Partition<NodeID> vp(g.num_nodes());                          // Divide nodes between PEs
  pvector<ScoreT> scores(vp.max_width, init_score, true);       // Symmetric partitioned array of scores (one for each node)
  pvector<ScoreT> outgoing_contrib(vp.max_width, true);         // Symmetric partitioned array of contribs
  for (int iter = 0; iter < max_iters; iter++) {
    *error = 0;
    for (NodeID n = vp.start; n < vp.end; n++)                  // Divide processing of nodes
      outgoing_contrib[vp.local_pos(n)] = scores[vp.local_pos(n)] / g.out_degree(n);
    shmem_barrier_all();
    for (NodeID u = vp.start; u < vp.end; u++) {
      incoming_total = 0;
      for (NodeID v : g.in_neigh(u)) {
        shmem_getmem(&temp, outgoing_contrib.begin()+vp.local_pos(v), sizeof(ScoreT), vp.recv(v));  // temp = outgoing_contrib[v]
        incoming_total += temp;
      }
      old_score = scores[vp.local_pos(u)];
      scores[vp.local_pos(u)] = base_score + kDamp * incoming_total;
      *error += fabs(scores[vp.local_pos(u)] - old_score);
    }
    shmem_barrier_all();
    shmem_double_sum_to_all(error, error, 1, 0, 0, vp.npes, pWrk, pSync);      // Reduction : +
    if (*error < epsilon)
      break;
  }
  return scores;
}

void PrintTopScores(const Graph &g, const pvector<ScoreT> &scores) {
  vector<pair<NodeID, ScoreT>> score_pairs(g.num_nodes());
  for (NodeID n=0; n < g.num_nodes(); n++) {
    score_pairs[n] = make_pair(n, scores[n]);
  }
  int k = 5;
  vector<pair<ScoreT, NodeID>> top_k = TopK(score_pairs, k);
  k = min(k, static_cast<int>(top_k.size()));
  for (auto kvp : top_k)
    cout << kvp.second << ":" << kvp.first << endl;
}


// Prints result to file to be read by original verifier
bool PRVerifier(const Graph &g, const pvector<ScoreT> &scores,
                        double target_error) {
  Partition<NodeID> vp(g.num_nodes());
  int* PRINTER = (int *) shmem_malloc(sizeof(int));
  *PRINTER = 0;
  shmem_barrier_all();
  shmem_int_wait_until(PRINTER, SHMEM_CMP_EQ, vp.pe);           // wait until previous PE puts your pe # in PRINTER
  ofstream shmem_out;
  shmem_out.open("pr_output.txt", ios::app);
  for (NodeID n = vp.start; n < vp.end; n++) {
    shmem_out << std::setprecision (5) << scores[vp.local_pos(n)] << std::endl;
  }
  shmem_out.close();
  if (!(vp.pe == vp.npes-1))
    shmem_int_p(PRINTER, vp.pe+1, vp.pe+1);             // who's next?
  return true;
}


int main(int argc, char* argv[]) {
  CLPageRank cli(argc, argv, "pagerank", 1e-4, 20);
  if (!cli.ParseArgs())
    return -1;

  char size_env[] = "SMA_SYMMETRIC_SIZE=16G";
  putenv(size_env);

  shmem_init();

  static long pSync[SHMEM_REDUCE_SYNC_SIZE];
  static long pWrk[SHMEM_REDUCE_MIN_WRKDATA_SIZE];
  static double dbl_pWrk[SHMEM_REDUCE_MIN_WRKDATA_SIZE];

  for (int i = 0; i < SHMEM_REDUCE_SYNC_SIZE; i++)
    pSync[i] = SHMEM_SYNC_VALUE;
  for (int i = 0; i < SHMEM_REDUCE_MIN_WRKDATA_SIZE; i++) {
    pWrk[i] = SHMEM_SYNC_VALUE;
    dbl_pWrk[i] = SHMEM_SYNC_VALUE;
  }

  static long* PRINT_LOCK = (long *) shmem_calloc(1, sizeof(long));

  {
    Builder b(cli, cli.do_verify());
    Graph g = b.MakeGraph(pWrk, pSync);
    auto PRBound = [&cli] (const Graph &g) {
      return PageRankPull(g, cli.max_iters(), pSync, dbl_pWrk, cli.tolerance());
    };
    auto VerifierBound = [&cli] (const Graph &g, const pvector<ScoreT> &scores) {
      return PRVerifier(g, scores, cli.tolerance());
    };
    BenchmarkKernel(cli, g, PRBound, PrintTopScores, VerifierBound);
  }

  shmem_free(PRINT_LOCK);
  shmem_finalize();
  return 0;
}
