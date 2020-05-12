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


typedef EdgePair<short, short> Edge;
typedef pvector<Edge> EdgeList;
typedef SlidingQueue<Edge> EdgeQ;
typedef QueueBuffer<Edge> QB;
typedef BuilderBase<int, int, int> Bob;

/*pvector<NodeID> InitParent(const Graph &g) {
  pvector<NodeID> parent(g.num_nodes(), true);
  #pragma omp parallel for
  for (NodeID n=0; n < g.num_nodes(); n++)
    parent[n] = g.out_degree(n) != 0 ? -g.out_degree(n) : -1;
  return parent;
}*/

// Partition parent array ~evenly across PEs
// Accessing node v on PE p means accessing node (n/k)*p + v in a complete parent array of n nodes and k PEs
pvector<NodeID> InitParent(const Graph &g) {
  int start, end;
  int pe = shmem_my_pe();
  int npes = shmem_n_pes();
  int offset = g.num_nodes()/npes;
  start = offset * pe;
  if (pe == npes-1) {
    end = g.num_nodes();  
  } else {
    end = start + offset; 
  }
  printf("N nodes: %ld | PE: %d | Start: %d | End: %d\n", g.num_nodes(), pe, start, end);
  pvector<NodeID> parent((end-start), true);
  #pragma omp parallel for
  for (NodeID n=start; n < end; n++)
    parent[n-start] = g.out_degree(n) != 0 ? -g.out_degree(n) : -1;
  return parent;
}



int main(int argc, char* argv[]) {
  CLApp cli(argc, argv, "breadth-first search");
  if (!cli.ParseArgs())
      return -1;

  static long FRONTIER_LOCK = 0;                                                      // Create a mutex lock in symmetric memory to control access to the frontier
  static long long pwrk[SHMEM_REDUCE_MIN_WRKDATA_SIZE];
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

  {
    int npes = shmem_n_pes();
    int pe = shmem_my_pe();
    Bitmap b(8, true);
    /*if (pe == 0) {
      b.set_bit(1);
      b.set_bit(2);
      b.set_bit(2);
    } else {
      b.set_bit(6);
    }*/
    b.reset();
    shmem_barrier_all();
    for (int i = 0; i < 8; i++)
      printf("PE: %d | Bit %d: %d\n", pe, i, b.get_bit(i));
    shmem_barrier_all();
    printf("Post-merge\n");
    b.merge(pwrk, pSync);
    shmem_barrier_all();
    for (int i = 0; i < 8; i++)
      printf("PE: %d | Bit %d: %d\n", pe, i, b.get_bit(i));
  }
/*  {
    void* builder_alloc = shmem_malloc(sizeof(Bob));
    Bob* builder = new(builder_alloc) Bob{cli};
    CSRGraph<int, int> g = builder->MakeGraph();
    //g.PrintTopology();
    shmem_barrier_all();
    pvector<NodeID> parent = InitParent(g);
    int i = 0;
    for (auto it = parent.begin(); it < parent.end(); it++){
      i++;
      printf("PE: %d | (%p => %d)\n", pe, (void *) it, *it);
    }
    shmem_free(builder);
  }*/                                                                                            // Extra scope to trigger deletion of graph, otherwise shmem destructor is screwy
  shmem_finalize();
  return 0;
}

void frontierior(void)
{
  int x = 5;
    /*void* frontier_alloc = shmem_malloc(sizeof(EdgeQ));                                // Just to emphasize that malloc only creates space, placement new concept in c++
    EdgeQ* frontier = new(frontier_alloc) EdgeQ{40};                                   // whereas new(frontier_buffer) EdgeQ is needed to call the edgeQ constructor
    size_t sighs = 3;
    QB* qb = new QB(*frontier, pe, npes, &FRONTIER_LOCK, sighs);
    for (auto it = EL->begin(); it < EL->end(); it++) {
        qb->push_back((*it));
    }
    for (int i = 0; i < sighs; i++){
        qb->push_back(Edge(pe+10, pe+10));
    }
    shmem_barrier_all();
    qb->flush();
    shmem_barrier_all();
    frontier->slide_window();
    shmem_barrier_all();
    frontier->~EdgeQ;
    */

}






 
