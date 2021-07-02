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

using namespace std;

int main(int argc, char* argv[]) {
  CLApp cli(argc, argv, "breadth-first search");
  if (!cli.ParseArgs())
    return -1;

  shmem_init();

  long* PRINT_LOCK = (long*) shmem_calloc(1, sizeof(long));

  static long pSync[SHMEM_REDUCE_SYNC_SIZE];
  static long pWrk[SHMEM_REDUCE_MIN_WRKDATA_SIZE];

  for (int i = 0; i < SHMEM_REDUCE_SYNC_SIZE; i++)
    pSync[i] = SHMEM_SYNC_VALUE;
  for (int i = 0; i < SHMEM_REDUCE_MIN_WRKDATA_SIZE; i++)
    pWrk[i] = SHMEM_SYNC_VALUE;

//  char size_env[] = "SHMEM_SYMMETRIC_SIZE=2048M";
//  putenv(size_env);

  {
    Builder b(cli);
    Graph g = b.MakeGraph(pWrk, pSync);
    g.PrintTopology();
  }

  shmem_free(PRINT_LOCK);
  shmem_finalize();
  return 0;
}
