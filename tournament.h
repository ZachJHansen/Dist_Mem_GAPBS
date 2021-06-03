  #ifndef TOURNAMENT_H_
  #define TOURNAMENT_H_

  #include <limits>
  #include "pvector.h"
  
  #define Parent(i) (i>>1)
  #define Left(i) (i<<1)
  #define Right(i) ((i<<1)+1)


  class TournamentTree {
    typedef std::pair<int64_t, int> node_pe_pair;
    public:
      TournamentTree(int* initial) {
        pe = shmem_my_pe();
        npes = shmem_n_pes();
        float pe_exp = log(npes) / log(2);             // log base 2 of npes
        if (floor(pe_exp) == ceil(pe_exp)) {            // npes is a power of 2
          complete = true;
          leaves = npes;
        } else {
          int pe_exp_round = (int) ceil(pe_exp);
          complete = false;
          leaves = pow(2, pe_exp_round);
        } 
        int i = 1;
        nodes = 1;
        height = 0;
        while (i < leaves) {                            // calculate number of nodes in complete tree
          nodes += leaves/i;
          i = i * 2;
          height++;
        } 
        tree = (node_pe_pair *) shmem_calloc(nodes+1, sizeof(node_pe_pair)); // Allocate enough space for complete tree with 1-indexing for convenience
        for (int i = 0; i < leaves-npes; i++)                                         // append infinite filler values to end of tree
          tree[nodes-i] = std::make_pair(std::numeric_limits<int64_t>::max(), NULL);
        for (int i = npes-1; i >= 0; i--)                                               // fill in remaining leaves
          tree[nodes-npes+i] = std::make_pair(initial[i], i);
        int64_t index = nodes;
        for (int h = height; h >= 0; h--) {
          for (int64_t i = 0; i < pow(2, h); i++)
            tree[Parent(index-i)] = min(tree[index-i+1], tree[index-i], std::less<node_pe_pair>());
          index -= pow(2, h);
        }

        /*for (int i = 0; i < leaves-npes; i++) {                                         // append infinite filler values to end of tree
          tree[nodes-i] = std::make_pair(std::numeric_limits<int64_t>::max(), NULL);
          i++;
          tree[nodes-i] = std::make_pair(std::numeric_limits<int64_t>::max(), NULL);
          tree[Parent(nodes-i)] = std::make_pair(std::numeric_limits<int64_t>::max(), NULL);
        }
        for (int i = leaves-npes; i < npes; i++) {                                              // fill in remaining leaves
          tree[npes-i] = std::make_pair(initial[i], i);
          i++;
          tree[npes-i] = std::make_pair(initial[i], i);
          tree[Parent(npes-i)] = min(tree[i-1], tree[i], std::greater<node_pe_pair>());
        }*/
      }

    void PrintTree() {
      for (int i = 1; i <= nodes; i++)
        printf("Value: %lu | PE: %d\n", tree[i].first, tree[i].second);
    }

    private:
      int pe;
      int npes; 
      int nodes;
      int leaves;
      int height;
      node_pe_pair* tree;
      bool complete;
  };

  #endif  // TOURNAMENT_H_
