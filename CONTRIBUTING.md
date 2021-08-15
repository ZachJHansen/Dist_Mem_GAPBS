We are looking for contributors! This is a big project, with many areas to improve upon. Areas to be improved or experimented with include:

Graph Partitioning - One limitation of the current vertex partitioning approach is that graphs with skewed degree distributions waste a significant amount of symmetric heap space. Alternative partitioning schemes should be explored.

Alternative Kernels - Our prototype implementation is functional, but not optimized. We hope to add implementations of alternative algorithms that are specifically designed for distributed memory settings.

Heuristic Tuning - The performance of SSSP, BFS, and TC rely heavily on finding effective parameters for their heuristics. We have not yet explored how adapting the GAPBS to a distributed memory model affects these heuristics. 

If you have any questions/suggestions/issues, please open an issue or email (zachhansen@unomaha.edu).
