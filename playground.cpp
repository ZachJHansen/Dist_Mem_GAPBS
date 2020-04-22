#include "bitmap.h"
#include "builder.h"
#include "command_line.h"
#include "graph.h"
#include "platform_atomics.h"
#include "pvector.h"
#include "sliding_queue.h"
#include "timer.h"

typedef EdgePair<int, int> Edge;
typedef pvector<Edge> EdgeList;
typedef SlidingQueue<Edge> EdgeQ;
typedef QueueBuffer<Edge> QB;
typedef BuilderBase<int, int, int> Bob;

int main(int argc, char* argv[]) {
  CLApp cli(argc, argv, "breadth-first search");
  if (!cli.ParseArgs())
      return -1;

  static long FRONTIER_LOCK = 0;                                                      // Create a mutex lock in symmetric memory to control access to the frontier
  static long PRINT_LOCK = 0;
  static long pSync[SHMEM_COLLECT_SYNC_SIZE];
  for (int i = 0; i < SHMEM_COLLECT_SYNC_SIZE; i++)
    pSync[i] = SHMEM_SYNC_VALUE;    

  shmem_init();

  int npes = shmem_n_pes();
  int pe = shmem_my_pe();
  int scale = 3;
  int degree = 2;

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
    */

  {
    void* builder_alloc = shmem_malloc(sizeof(Bob));
    Bob* builder = new(builder_alloc) Bob{cli};
    CSRGraph<int, int> g = builder->MakeGraph();
    g.PrintTopology();
    shmem_barrier_all();
    shmem_free(builder);
  }                                                                                             // Graph must go out of scope before shmem_finalize for destructor to work properly


// Now we need a parent array in symmetric memory
//    frontier->~EdgeQ();
 //   shfree(frontier_alloc);
  shmem_finalize();
  return 0;
}



