  #ifndef TOURNAMENT_H_
  #define TOURNAMENT_H_

  #include <limits>
  #include <tuple>
  #include "pvector.h"
  
  #define Parent(i) (i>>1)
  #define Left(i) (i<<1)
  #define Right(i) ((i<<1)+1)

  class TournamentTree {
    public:
      typedef std::tuple<int64_t, int, int> deg_node_pe;          // <Degree, NodeID_, PE #>
      typedef std::pair<int64_t, int> degree_node_p;              // replace int with NodeID_ eventually
      std::vector<int> list_lengths;

      TournamentTree() {}

      TournamentTree(degree_node_p* initial, const pvector<degree_node_p> &partially_sorted_list) {
        pe = shmem_my_pe();
        npes = shmem_n_pes();
        partial = partially_sorted_list.begin();
        pop_counter = (int *) shmem_malloc(sizeof(int));
        *pop_counter = 1;                               // Every PE contributes the first element of their partial list
        float pe_exp = log(npes) / log(2);              // log base 2 of npes
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
        empty_leaf = (degree_node_p *) shmem_malloc(sizeof(degree_node_p)); 
        tree = (deg_node_pe *) shmem_calloc(nodes+1, sizeof(deg_node_pe));              // Allocate enough space for complete tree with 1-indexing for convenience
        //printf("Leaves: %d | Nodes: %d | Height: %d\n", leaves, nodes, height);
        for (int i = 0; i < leaves-npes; i++)                                           // append infinite filler values to end of tree
          tree[nodes-i] = std::make_tuple(std::numeric_limits<int64_t>::max(), 0, 0);
        for (int i = 0; i < npes; i++) {                                              // fill in remaining leaves, record list lengths
          tree[nodes-leaves+1+i] = std::make_tuple(initial[i].first, initial[i].second, i);
          list_lengths.push_back(partially_sorted_list.local_width(i));
        }
        int64_t index = nodes;                                                          // fill tree from leaves to root
        for (int h = height; h >= 0; h--) {
          for (int64_t i = 0; i < pow(2, h); i++)
            tree[Parent(index-i)] = min(tree[index-i+1], tree[index-i], std::less<deg_node_pe>());
          index -= pow(2, h);
        }
        printf("PE %d | tree: %p\n", pe, (void*) tree);
        shmem_barrier_all();
      }

    void print_tree() {
      for (int i = 1; i <= nodes; i++)
        printf("Value: %lu | PE: %d\n", std::get<0>(tree[i]), std::get<1>(tree[i]));
    }

    // Transfer the contents of the tree from this PE to the next PE
    void transfer(int current_pe) {
      printf("PE %d | tree: %p\n", pe, (void*) tree);
      size_t bytes = (nodes+1) * sizeof(deg_node_pe);
      shmem_putmem(tree, tree, bytes, current_pe+1);
      shmem_quiet();
    }

    degree_node_p pop_root() {
      int64_t val = std::get<0>(tree[1]);
      int node = std::get<1>(tree[1]);
      int loser = std::get<2>(tree[1]);                                                       // the pe that supplied the root node 
      int64_t index = (int64_t) shmem_int_atomic_fetch_inc(pop_counter, loser);                 // keep track of how deep in loser's list you are 
      int64_t node_id = nodes-leaves+1+loser;
      if (*pop_counter+1 < list_lengths[loser]) {                                             // Check to ensure we haven't emptied the list of the losing PE
        shmem_getmem(empty_leaf, partial+(index+1), sizeof(degree_node_p), loser);                          // request the next value from the pe that supplied the root node
        tree[node_id] = std::make_tuple(empty_leaf->first, empty_leaf->second, loser);
      } else {
        tree[node_id] = std::make_tuple(std::numeric_limits<int64_t>::max(), 0, loser);
      } 
      for (int h = height; h >= 0; h--) {                                               // Rebuild missing parts of the tree
        if (node_id % 2 == 0) {                                                                         // Even implies left child
          tree[Parent(node_id)] = min(tree[node_id], tree[node_id+1], std::less<deg_node_pe>());       // Compare with right sibling
        } else {
          tree[Parent(node_id)] = min(tree[node_id], tree[node_id-1], std::less<deg_node_pe>());       // Compare with left sibling
        }
        node_id = Parent(node_id);                                                      // Crawl up tree
      }
      return(std::make_pair(val, node));
    }

    private:
      int pe;
      int npes; 
      int nodes;
      int leaves;
      int height;
      int* pop_counter;
      deg_node_pe* tree;
      degree_node_p* empty_leaf;
      pvector<degree_node_p>::iterator partial; 
      bool complete;
  };

  #endif  // TOURNAMENT_H_
