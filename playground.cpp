#include "/home/zach/projects/GAPBS/gapbs/src/generator.h"
#include "sliding_queue.h"
#include <utility>

#define True 1
#define False 0

typedef EdgePair<int, int> Edge;
typedef pvector<Edge> EdgeList;
typedef SlidingQueue<Edge> EdgeQ;
typedef QueueBuffer<Edge> QB;

int main(void)
{
    static long FRONTIER_LOCK = 0;                                                      // Create a mutex lock in symmetric memory to control access to the frontier
    shmem_init();

    int npes = shmem_n_pes();
    int pe = shmem_my_pe();
    int scale = 3;
    int degree = 2;

    int *num_nodes = (int*) shmem_malloc(sizeof(int));
    size_t *num_edges = (size_t*) shmem_malloc(sizeof(size_t));
    *num_nodes = 1l << scale;                                                           // N = 2^scale nodes
    *num_edges = *num_nodes * degree;
    EdgeList* EL = (EdgeList*) shmem_calloc(*num_edges, sizeof(Edge));                  // Create space for an edge list in symmetric memory across all PEs.

    long pSync[_SHMEM_BCAST_SYNC_SIZE];                                                 // Initialize symmetric work array pSync
    size_t *el_size = (size_t*) shmem_malloc(sizeof(size_t));                           // Create a size_t slot in symmetric memory called el_size

    Generator<int, int, int> genuine{scale, degree};                                    // Generate an unweighted graph with 2^scale vertices, each with degree connections?
    *EL = genuine.GenerateEL(True);                                                     // Generate a uniform edge list in symmetric memory

    void* frontier_alloc = shmem_malloc(sizeof(EdgeQ));                                // Just to emphasize that malloc only creates space, placement new concept in c++
    EdgeQ* frontier = new(frontier_alloc) EdgeQ{40};                                   // whereas new(frontier_buffer) EdgeQ is needed to call the edgeQ constructor
/*    if (pe == 0){
        for (auto it = (*EL).begin(); it < (*EL).end(); it++) {
            frontier->push_back((*it));
        }
        frontier->slide_window();                                                       // foo->method() is the same as (*foo).method!
        for (auto it = frontier->begin(); it < frontier->end(); it++) {
            printf("PE: %d | Src: %d | Dest: %d\n", pe, (*it).u, (*it).v);
        }
    }*/
    //EdgeQ& f = *frontier;
    size_t sighs = 3;
    QB* qb = new QB(*frontier, pe, npes, &FRONTIER_LOCK, sighs);
    /*for (auto it = EL->begin(); it < EL->end(); it++) {
        qb->push_back((*it));
    }*/
    for (int i = 0; i < sighs; i++){
        qb->push_back(Edge(pe+10, pe+10));
    }
    shmem_barrier_all();
    qb->flush();
    shmem_barrier_all();
    frontier->slide_window();
    shmem_barrier_all();
    if (pe == 2){
        for (auto it = frontier->begin(); it < frontier->end(); it++) {
            printf("PE: %d | Src: %d | Dest: %d\n", pe, (*it).u, (*it).v);
        }
    }
    shmem_barrier_all();
    if (pe == 4){
        for (auto it = frontier->begin(); it < frontier->end(); it++) {
            printf("PE: %d | Src: %d | Dest: %d\n", pe, (*it).u, (*it).v);
        }
    }
    //qb->push_back(*(EL->begin()));
    
    // Next step is shmemifying the queue buffers (currently each thread manages a queuebuffer, we need each PE to manage a qb) 
    // Now we need a parent array in symmetric memory
    frontier->~EdgeQ();
    shfree(frontier_alloc);
    shfree(num_nodes);
    shfree(num_edges);
    shfree(EL);                                                         // Deallocate the edge_list on all PEs
    shfree(el_size);
    shmem_finalize();
    return 0;
}
