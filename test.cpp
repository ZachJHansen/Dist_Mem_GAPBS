#include "shmem.h"
#include "pvector.h"
#include "sliding_queue.h"


void DOBFS(long *FL)
{
  SlidingQueue<int> parent(10);
  QueueBuffer<int> buff(parent, FL);
}

int main(void)
{
  static long FL = 0;
  static int pWrk[SHMEM_REDUCE_MIN_WRKDATA_SIZE];
  static long pSync[SHMEM_REDUCE_SYNC_SIZE];
  for (int i = 0; i < SHMEM_REDUCE_SYNC_SIZE; i++)
    pSync[i] = SHMEM_SYNC_VALUE;
  for (int i = 0; i < SHMEM_REDUCE_MIN_WRKDATA_SIZE; i++)
    pWrk[i] = SHMEM_SYNC_VALUE;
 
  shmem_init();

  int pe, npes;
  pe = shmem_my_pe();
  npes = shmem_n_pes();
  int* scout_count = (int *) shmem_malloc(sizeof(int));
  int* local_scout = (int *) shmem_malloc(sizeof(int));
//  DOBFS(&FL);
  *local_scout = pe;
  int master[10] = {0,1,2,3,4,5,6,7,8,9};
  for (i : master) {
    int p = i/3;
    if (p > 2) {
      printf("special i: %d | PE: %d | Local: %d\n", i, 2, i - ((10*2)/3));
    } else {
      printf("i: %d | PE: %d | Local: %d\n", i, p, i%3);
    }
  }
//  printf("PE: %d - Scout (local): %d\n", pe, *local_scout);
  //shmem_int_sum_to_all(scout_count, local_scout, 1, 0, 0, npes, pWrk, pSync);
  //printf("PE: %d - Scout (total): %d\n", pe, *scout_count);
  
  shmem_finalize();

  return 0;
}
