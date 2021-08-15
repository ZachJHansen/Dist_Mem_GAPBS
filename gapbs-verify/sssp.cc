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
Author: Scott Beamer, Yunming Zhang

Returns array of distances for all vertices from given source vertex

This SSSP implementation makes use of the ∆-stepping algorithm [1]. The type
used for weights and distances (WeightT) is typedefined in benchmark.h. The
delta parameter (-d) should be set for each input graph. This implementation
incorporates a new bucket fusion optimization [2] that significantly reduces
the number of iterations (& barriers) needed.

The bins of width delta are actually all thread-local and of type std::vector
so they can grow but are otherwise capacity-proportional. Each iteration is
done in two phases separated by barriers. In the first phase, the current
shared bin is processed by all threads. As they find vertices whose distance
they are able to improve, they add them to their thread-local bins. During this
phase, each thread also votes on what the next bin should be (smallest
non-empty bin). In the next phase, each thread copies its selected
thread-local bin into the shared bin.

Once a vertex is added to a bin, it is not removed, even if its distance is
later updated and it now appears in a lower bin. We find ignoring vertices if
their distance is less than the min distance for the current bin removes
enough redundant work to be faster than removing the vertex from older bins.

The bucket fusion optimization [2] executes the next thread-local bin in
the same iteration if the vertices in the next thread-local bin have the
same priority as those in the current shared bin. This optimization greatly
reduces the number of iterations needed without violating the priority-based
execution order, leading to significant speedup on large diameter road networks.

[1] Ulrich Meyer and Peter Sanders. "δ-stepping: a parallelizable shortest path
    algorithm." Journal of Algorithms, 49(1):114–152, 2003.

[2] Yunming Zhang, Ajay Brahmakshatriya, Xinyi Chen, Laxman Dhulipala,
    Shoaib Kamil, Saman Amarasinghe, and Julian Shun. "Optimizing ordered graph
    algorithms with GraphIt." The 18th International Symposium on Code Generation
    and Optimization (CGO), pages 158-170, 2020.
*/


using namespace std;

const WeightT kDistInf = numeric_limits<WeightT>::max()/2;
const size_t kMaxBin = numeric_limits<size_t>::max()/2;
const size_t kBinSizeThreshold = 1000;

inline
void RelaxEdges(const WGraph &g, NodeID u, WeightT delta,
                pvector<WeightT> &dist, vector <vector<NodeID>> &local_bins) {
  //printf("Relaxing edge %d\n", u);
  for (WNode wn : g.out_neigh(u)) {
    WeightT old_dist = dist[wn.v];
    WeightT new_dist = dist[u] + wn.w;
    //printf("U: %d | old_dist: %d | new_dist: %d\n", u, old_dist, new_dist);
    while (new_dist < old_dist) {
      if (compare_and_swap(dist[wn.v], old_dist, new_dist)) {
        size_t dest_bin = new_dist/delta;
        if (dest_bin >= local_bins.size())
          local_bins.resize(dest_bin+1);
        local_bins[dest_bin].push_back(wn.v);
        break;
      }
      old_dist = dist[wn.v];      // swap failed, recheck dist update & retry
    }
  }
}

pvector<WeightT> DeltaStep(const WGraph &g, NodeID source, WeightT delta) {
  printf("Source: %d\n", source);
  Timer t;
  pvector<WeightT> dist(g.num_nodes(), kDistInf);
  dist[source] = 0;
  pvector<NodeID> frontier(g.num_edges_directed());
  // two element arrays for double buffering curr=iter&1, next=(iter+1)&1
  size_t shared_indexes[2] = {0, kMaxBin};
  size_t frontier_tails[2] = {1, 0};
  frontier[0] = source;
  t.Start();
  #pragma omp parallel
  {
    vector<vector<NodeID> > local_bins(0);
    size_t iter = 0;
    while (shared_indexes[iter&1] != kMaxBin) {
      size_t &curr_bin_index = shared_indexes[iter&1];
      size_t &next_bin_index = shared_indexes[(iter+1)&1];
      size_t &curr_frontier_tail = frontier_tails[iter&1];
      size_t &next_frontier_tail = frontier_tails[(iter+1)&1];
      //printf("While start | iter = %lu, shared_index[iter&1] = %lu\n", iter, shared_indexes[(iter)&1]);
      //printf("curr_bin_index: %lu | next_bin_index: %lu | curr_front_tail: %lu | next_front_tail: %lu \n", curr_bin_index, next_bin_index, curr_frontier_tail, next_frontier_tail);
      //printf("iter = %lu, iter&1 = %lu, (iter+1)&1 = %lu\n", iter, iter&1, (iter+1)&1);
      #pragma omp for nowait schedule(dynamic, 64)
      //printf("Current frontier tail is : %lu\n", curr_frontier_tail);
      for (size_t i=0; i < curr_frontier_tail; i++) {
        NodeID u = frontier[i];
        //printf("U = %d\n", u);

        //printf("i: %d | Udist: %d\n", i, dist[u]);
        //printf("dist[%d] = %d, curr_bin_index = %d, delta = %d\n", u, dist[u], static_cast<WeightT>(curr_bin_index), delta);
        if (dist[u] >= delta * static_cast<WeightT>(curr_bin_index)) {
          //printf("Relaxing opp 1\n");
          RelaxEdges(g, u, delta, dist, local_bins);
        }
      }
      while (curr_bin_index < local_bins.size() &&
             !local_bins[curr_bin_index].empty() &&
             local_bins[curr_bin_index].size() < kBinSizeThreshold) {
        //printf("mystery\n");
        vector<NodeID> curr_bin_copy = local_bins[curr_bin_index];
        local_bins[curr_bin_index].resize(0);
        //printf("Post resize relaxation\n");
        for (NodeID u : curr_bin_copy) {
          //printf("Relaxing opp 2: u = %d\n", u);
          RelaxEdges(g, u, delta, dist, local_bins);
        }
      }
      //if (curr_bin_index < local_bins.size())
        //printf("check 1\n");
      for (size_t i=curr_bin_index; i < local_bins.size(); i++) {
        if (!local_bins[i].empty()) {
          #pragma omp critical
          next_bin_index = min(next_bin_index, i);
          //printf("next: %d\n", next_bin_index);
          break;
        }
      }
      #pragma omp barrier
      //printf("The chosen next bin is %lu\n", next_bin_index);
      #pragma omp single nowait
      {
        //printf("PE %d | next chosen bin is %lu\n", 0, shared_indexes[((iter)+1)&1]);
        t.Stop();
        PrintStep(curr_bin_index, t.Millisecs(), curr_frontier_tail);
        t.Start();
        curr_bin_index = kMaxBin;
        curr_frontier_tail = 0;
      }
      //printf("Curr bin index = %lu | curr_front = %lu | shared_indexes[iter&1] = %lu | shared_indexes[(iter+1)&1] = %lu\n", curr_bin_index, curr_frontier_tail, shared_indexes[iter&1], shared_indexes[(iter+1)&1]);


      //printf("Prior to iter inc | iter = %lu, shared_index[iter&1] = %lu, curr_front = %lu\n", iter, shared_indexes[(iter)&1], curr_frontier_tail);
      if (next_bin_index < local_bins.size()) {
  /*      printf("Frontier contents b4 copy: ");
        for (WeightT t : frontier)
          printf("%d ", t);
        printf("\n");*/
        size_t copy_start = fetch_and_add(next_frontier_tail,
                                          local_bins[next_bin_index].size());
        ofstream out;
        out.open("/home/zach/projects/Dist_Mem_GAPBS/Dist_Mem_GAPBS/counts.txt", ios::app);
        out << iter << " " << local_bins[next_bin_index].size() << endl;
        out.close();
        //printf("end - begin = %d\n", local_bins[next_bin_index].end() - local_bins[next_bin_index].begin());
        copy(local_bins[next_bin_index].begin(),
             local_bins[next_bin_index].end(), frontier.data() + copy_start);
        local_bins[next_bin_index].resize(0);
        /*printf("Frontier contents after copy: ");
        for (WeightT t : frontier)
          printf("%d ", t);
        printf("\n");*/
      }
      iter++;

      //printf("Prior to iter inc | iter = %lu, shared_index[iter&1] = %lu, curr_front = %lu\n", iter, shared_indexes[(iter)&1], curr_frontier_tail);
      #pragma omp barrier
      //printf("While end | iter = %lu, shared_index[iter&1] = %lu\n", iter, shared_indexes[(iter)&1]);
    }
    #pragma omp single
    cout << "took " << iter << " iterations" << endl;
  }
  //for (int i = 0; i < dist.size(); i++)
    //printf("dist[%d] = %d\n", i, dist[i]);
  return dist;
}


void PrintSSSPStats(const WGraph &g, const pvector<long> &dist) {
  auto NotInf = [](WeightT d) { return d != kDistInf; };
  int64_t num_reached = count_if(dist.begin(), dist.end(), NotInf);
  cout << "SSSP Tree reaches " << num_reached << " nodes" << endl;
}


// Compares against simple serial implementation
bool SSSPVerifier(const WGraph &g, NodeID source,
                  const pvector<long> &dist_to_test) {
  // Serial Dijkstra implementation to get oracle distances

  printf("Source: %d\n", source);
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
  WeightedBuilder b(cli);
  WGraph g = b.MakeGraph();
  g.PrintTopology(true);
  g.PrintTopology(false);
  SourcePicker<WGraph> sp(g, cli.start_vertex());
  auto SSSPBound = [&sp, &cli] (const WGraph &g) {
    return DeltaStep(g, sp.PickNext(), cli.delta());
  };
  /*pvector<WeightT> parent = SSSPBound(g);
  int i = 0;
  for (WeightT p : parent)
    printf("dist(%d) = %lu\n", i++, p);*/
  SourcePicker<WGraph> vsp(g, cli.start_vertex());
  //auto discard = vsp.PickNext();
  auto VerifierBound = [&vsp] (const WGraph &g, const pvector<long> &dist) {
    return SSSPVerifier(g, vsp.PickNext(), dist);
  };
  //pvector<WeightT> dist = SSSPBound(g);
  SSSPBenchmarkKernel(cli, g, SSSPBound, PrintSSSPStats, VerifierBound);
  return 0;
}
