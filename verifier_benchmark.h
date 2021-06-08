// Copyright (c) 2015, The Regents of the University of California (Regents)
// See LICENSE.txt for license details

#ifndef BENCHMARK_H_
#define BENCHMARK_H_

#include <algorithm>
#include <cinttypes>
#include <functional>
#include <random>
#include <utility>
#include <vector>

#include "builder.h"
#include "graph.h"
#include "timer.h"
#include "util.h"
#include "writer.h"


/*
GAP Benchmark Suite
File:   Benchmark
Author: Scott Beamer

Various helper functions to ease writing of kernels
*/


// Default type signatures for commonly used types
typedef int32_t NodeID;
typedef int32_t WeightT;
typedef NodeWeight<NodeID, WeightT> WNode;

typedef CSRGraph<NodeID> Graph;
typedef CSRGraph<NodeID, WNode> WGraph;

typedef BuilderBase<NodeID, NodeID, WeightT> Builder;
typedef BuilderBase<NodeID, WNode, WeightT> WeightedBuilder;

typedef WriterBase<NodeID, NodeID> Writer;
typedef WriterBase<NodeID, WNode> WeightedWriter;


// Used to pick random non-zero degree starting points for search algorithms
template<typename GraphT_>
class SourcePicker {
 public:
  explicit SourcePicker(const GraphT_ &g, NodeID given_source = -1)
      : given_source(given_source), rng(kRandSeed), udist(0, g.num_nodes()-1),
        g_(g) {}

  NodeID PickNext() {
    if (given_source != -1)
      return given_source;
    NodeID source;
    do {
      source = udist(rng);
    } while (g_.out_degree(source) == 0);
    return source;
  }

 private:
  NodeID given_source;
  std::mt19937 rng;
  std::uniform_int_distribution<NodeID> udist;
  const GraphT_ &g_;
};


// Returns k pairs with largest values from list of key-value pairs
template<typename KeyT, typename ValT>
std::vector<std::pair<ValT, KeyT>> TopK(
    const std::vector<std::pair<KeyT, ValT>> &to_sort, size_t k) {
  std::vector<std::pair<ValT, KeyT>> top_k;
  ValT min_so_far = 0;
  for (auto kvp : to_sort) {
    if ((top_k.size() < k) || (kvp.second > min_so_far)) {
      top_k.push_back(std::make_pair(kvp.second, kvp.first));
      std::sort(top_k.begin(), top_k.end(),
                std::greater<std::pair<ValT, KeyT>>());
      if (top_k.size() > k)
        top_k.resize(k);
      min_so_far = top_k.back().first;
    }
  }
  return top_k;
}


bool VerifyUnimplemented(...) {
  std::cout << "** verify unimplemented **" << std::endl;
  return false;
}

// Verify that the distributed memory benchmark produced the expected result (recorded in pr_output.txt)
// by building and processing the same graph with the original implementation and comparing with pr_output.txt
// No timing of this secondary build-process-verify procedure is recorded 
template<typename GraphT_, typename GraphFunc, typename AnalysisFunc,
         typename VerifierFunc>
void PRBenchmarkKernel(const CLApp &cli, const GraphT_ &g,
                     GraphFunc kernel, AnalysisFunc stats,
                     VerifierFunc verify) {
  pvector<float> master(g.num_nodes() * cli.num_trials());
  master.clear();
  std::ifstream shmem_out;
  shmem_out.open("/home/zhansen/Dist_Mem_GAPBS/pr_output.txt"); //open a file to perform read operation using file object
  if (shmem_out.is_open()){   //checking whether the file is open
    std::string tp;
    while(getline(shmem_out, tp))  //read data from file object and put it into string.
        master.push_back(std::stof(tp));        // page rank produces an array of scores (floats)
    shmem_out.close();   //close the file object.
  }  
  std::remove("/home/zhansen/Dist_Mem_GAPBS/pr_output.txt");               // delete file
  g.PrintStats();
  pvector<float> result(g.num_nodes());
  for (int iter=0; iter < cli.num_trials(); iter++) {
    result.clear();
    pvector<float> original_scores = kernel(g);
    for (int i = iter*g.num_nodes(); i < (iter+1)*g.num_nodes(); i++)
      result.push_back(master[i]);
    if (cli.do_analysis() && (iter == (cli.num_trials()-1)))
      stats(g, result);
    if (cli.do_verify()) {
      PrintLabel("Verification",
                 verify(std::ref(original_scores), std::ref(result)) ? "PASS" : "FAIL");
    }
  }
}


template<typename GraphT_, typename GraphFunc, typename AnalysisFunc,
         typename VerifierFunc>
void BFSBenchmarkKernel(const CLApp &cli, const GraphT_ &g,
                     GraphFunc kernel, AnalysisFunc stats,
                     VerifierFunc verify) {
  printf("These times aren't valid, I'm just checking the shmem_output.txt for validity\n");
  pvector<int> master(g.num_nodes() * cli.num_trials());
  master.clear();
  std::ifstream shmem_out;
  shmem_out.open("/home/zhansen/Dist_Mem_GAPBS/bfs_output.txt"); //open a file to perform read operation using file object
  if (shmem_out.is_open()){                     //checking whether the file is open
    std::string tp;
    while(getline(shmem_out, tp))               //read data from file object and put it into string.
        master.push_back(std::stoi(tp));        // bfs and connected components produce parent arrays (ints)
    shmem_out.close();                          //close the file object.
  } 
  std::remove("/home/zhansen/Dist_Mem_GAPBS/bfs_output.txt");               // delete file
  g.PrintStats();
  double total_seconds = 0;
  Timer trial_timer;
  pvector<int> result(g.num_nodes());
  for (int iter=0; iter < cli.num_trials(); iter++) {
    result.clear();
    trial_timer.Start();
    for (int i = iter*g.num_nodes(); i < (iter+1)*g.num_nodes(); i++)
      result.push_back(master[i]);
    
    trial_timer.Stop();
    PrintTime("Trial Time", trial_timer.Seconds());
    total_seconds += trial_timer.Seconds();
    if (cli.do_analysis() && (iter == (cli.num_trials()-1)))
      stats(g, result);
    if (cli.do_verify()) {
      trial_timer.Start();
      PrintLabel("Verification",
                 verify(std::ref(g), std::ref(result)) ? "PASS" : "FAIL");
      trial_timer.Stop();
      PrintTime("Verification Time", trial_timer.Seconds());
    }
  }
  PrintTime("Average Time", total_seconds / cli.num_trials());
}


template<typename GraphT_, typename GraphFunc, typename AnalysisFunc,
         typename VerifierFunc>
void CCBenchmarkKernel(const CLApp &cli, const GraphT_ &g,
                     GraphFunc kernel, AnalysisFunc stats,
                     VerifierFunc verify) {
  printf("These times aren't valid, I'm just checking the shmem_output.txt for validity\n");
  pvector<int> master(g.num_nodes() * cli.num_trials());
  master.clear();
  std::ifstream shmem_out;
  shmem_out.open("/home/zhansen/Dist_Mem_GAPBS/cc_output.txt"); //open a file to perform read operation using file object
  if (shmem_out.is_open()){                     //checking whether the file is open
    std::string tp;
    while(getline(shmem_out, tp))               //read data from file object and put it into string.
        master.push_back(std::stoi(tp));        // bfs and connected components produce parent arrays (ints)
    shmem_out.close();                          //close the file object.
  } 
  std::remove("/home/zhansen/Dist_Mem_GAPBS/cc_output.txt");               // delete file
  g.PrintStats();
  double total_seconds = 0;
  Timer trial_timer;
  pvector<int> result(g.num_nodes());
  for (int iter=0; iter < cli.num_trials(); iter++) {
    result.clear();
    trial_timer.Start();
    for (int i = iter*g.num_nodes(); i < (iter+1)*g.num_nodes(); i++)
      result.push_back(master[i]);
    
    trial_timer.Stop();
    PrintTime("Trial Time", trial_timer.Seconds());
    total_seconds += trial_timer.Seconds();
    if (cli.do_analysis() && (iter == (cli.num_trials()-1)))
      stats(g, result);
    if (cli.do_verify()) {
      trial_timer.Start();
      PrintLabel("Verification",
                 verify(std::ref(g), std::ref(result)) ? "PASS" : "FAIL");
      trial_timer.Stop();
      PrintTime("Verification Time", trial_timer.Seconds());
    }
  }
  PrintTime("Average Time", total_seconds / cli.num_trials());
}


template<typename GraphT_, typename GraphFunc, typename AnalysisFunc,
         typename VerifierFunc>
void TCBenchmarkKernel(const CLApp &cli, const GraphT_ &g,
                     GraphFunc kernel, AnalysisFunc stats,
                     VerifierFunc verify) {
  printf("These times aren't valid, I'm just checking the shmem_output.txt for validity\n");
  pvector<int> master(cli.num_trials());
  master.clear();
  std::ifstream shmem_out;
  shmem_out.open("/home/zhansen/Dist_Mem_GAPBS/tc_output.txt"); //open a file to perform read operation using file object
  if (shmem_out.is_open()){   //checking whether the file is open
    std::string tp;
    while(getline(shmem_out, tp)){  //read data from file object and put it into string.
      master.push_back(std::stoi(tp));
    }
    shmem_out.close();   //close the file object.
  } 
  std::remove("/home/zhansen/Dist_Mem_GAPBS/tc_output.txt");               // delete file
  g.PrintStats();
  double total_seconds = 0;
  Timer trial_timer;
  for (int iter=0; iter < cli.num_trials(); iter++) {
    size_t result = master[iter]; 
    trial_timer.Stop();
    PrintTime("Trial Time", trial_timer.Seconds());
    total_seconds += trial_timer.Seconds();
    if (cli.do_analysis() && (iter == (cli.num_trials()-1)))
      stats(g, result);
    if (cli.do_verify()) {
      trial_timer.Start();
      PrintLabel("Verification",
                 verify(std::ref(g), std::ref(result)) ? "PASS" : "FAIL");
      trial_timer.Stop();
      PrintTime("Verification Time", trial_timer.Seconds());
    }
  }
}

template<typename GraphT_, typename GraphFunc, typename AnalysisFunc,
         typename VerifierFunc>
void SSSPBenchmarkKernel(const CLApp &cli, const GraphT_ &g,
                     GraphFunc kernel, AnalysisFunc stats,
                     VerifierFunc verify) {
  printf("These times aren't valid, I'm just checking the shmem_output.txt for validity\n");
  pvector<long> master(g.num_nodes() * cli.num_trials());
  master.clear();
  std::ifstream shmem_out;
  shmem_out.open("/home/zhansen/Dist_Mem_GAPBS/sssp_output.txt"); //open a file to perform read operation using file object
  if (shmem_out.is_open()){                     //checking whether the file is open
    std::string tp;
    while(getline(shmem_out, tp))               //read data from file object and put it into string.
        master.push_back(std::stol(tp));        // bfs and connected components produce parent arrays (ints)
    shmem_out.close();                          //close the file object.
  }
  std::remove("/home/zhansen/Dist_Mem_GAPBS/sssp_output.txt");               // delete file
  g.PrintStats();
  double total_seconds = 0;
  Timer trial_timer;
  pvector<long> result(g.num_nodes());
  for (int iter=0; iter < cli.num_trials(); iter++) {
    result.clear();
    trial_timer.Start();
    for (int i = iter*g.num_nodes(); i < (iter+1)*g.num_nodes(); i++)
      result.push_back(master[i]);

    trial_timer.Stop();
    PrintTime("Trial Time", trial_timer.Seconds());
    total_seconds += trial_timer.Seconds();
    if (cli.do_analysis() && (iter == (cli.num_trials()-1)))
      stats(g, result);
    if (cli.do_verify()) {
      trial_timer.Start();
      PrintLabel("Verification",
                 verify(std::ref(g), std::ref(result)) ? "PASS" : "FAIL");
      trial_timer.Stop();
      PrintTime("Verification Time", trial_timer.Seconds());
    }
  }
  PrintTime("Average Time", total_seconds / cli.num_trials());
}

template<typename GraphT_, typename GraphFunc, typename AnalysisFunc,
         typename VerifierFunc>
void BCBenchmarkKernel(const CLApp &cli, const GraphT_ &g,
                     GraphFunc kernel, AnalysisFunc stats,
                     VerifierFunc verify) {
  printf("These times aren't valid, I'm just checking the shmem_output.txt for validity\n");
  pvector<float> master(g.num_nodes() * cli.num_trials());
  master.clear();
  std::ifstream shmem_out;
  shmem_out.open("/home/zhansen/Dist_Mem_GAPBS/bc_output.txt"); //open a file to perform read operation using file object
  if (shmem_out.is_open()){                     //checking whether the file is open
    std::string tp;
    while(getline(shmem_out, tp))               //read data from file object and put it into string.
        master.push_back(std::stof(tp));        // bfs and connected components produce parent arrays (ints)
    shmem_out.close();                          //close the file object.
  }
  std::remove("/home/zhansen/Dist_Mem_GAPBS/bc_output.txt");               // delete file
  g.PrintStats();
  double total_seconds = 0;
  Timer trial_timer;
  pvector<float> result(g.num_nodes());
  for (int iter=0; iter < cli.num_trials(); iter++) {
    result.clear();
    trial_timer.Start();
    for (int i = iter*g.num_nodes(); i < (iter+1)*g.num_nodes(); i++)
      result.push_back(master[i]);

    trial_timer.Stop();
    PrintTime("Trial Time", trial_timer.Seconds());
    total_seconds += trial_timer.Seconds();
    if (cli.do_analysis() && (iter == (cli.num_trials()-1)))
      stats(g, result);
    if (cli.do_verify()) {
      trial_timer.Start();
      PrintLabel("Verification",
                 verify(std::ref(g), std::ref(result)) ? "PASS" : "FAIL");
      trial_timer.Stop();
      PrintTime("Verification Time", trial_timer.Seconds());
    }
  }
  PrintTime("Average Time", total_seconds / cli.num_trials());
}



// Calls (and times) kernel according to command line arguments
template<typename GraphT_, typename GraphFunc, typename AnalysisFunc,
         typename VerifierFunc>
void BenchmarkKernel(const CLApp &cli, const GraphT_ &g,
                     GraphFunc kernel, AnalysisFunc stats,
                     VerifierFunc verify, char bench_type) {
  /*switch(bench_type) {
    case 'b':
      BFSBench(cli, g, kernel, stats, verify, 'b');
    case 't':
      TCBench(cli, g, kernel, stats, verify);
    case 's':
      SSSPBench(cli, g, kernel, stats, verify);
    case 'p':
      PRBench(cli, g, kernel, stats, verify, 'p');
    case 'c':
      BFSBench(cli, g, kernel, stats, verify, 'c');
    default:
      printf("Specify b, t, p, c, or s\n");
  }*/
  //TCBench(cli, g, kernel, stats, verify);
 // PRBench(cli, g, kernel, stats, verify, 'p');
   BFSBench(cli, g, kernel, stats, verify, 'b');
}


#endif  // BENCHMARK_H_
