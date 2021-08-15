// Copyright (c) 2015, The Regents of the University of California (Regents)
// See LICENSE.txt for license details

#ifndef UTIL_H_
#define UTIL_H_

#include <stdio.h>
#include <cinttypes>
#include <string>

#include "timer.h"


/*
GAP Benchmark Suite
Author: Scott Beamer

Miscellaneous helpers that don't fit into classes

Zach's Notes:
 - Only PE 0 prints time ops
*/


static const int64_t kRandSeed = 27491095;


void PrintLabel(const std::string &label, const std::string &val) {
  if (shmem_my_pe() == 0)
    printf("PE: %d | %-21s%7s\n", shmem_my_pe(), (label + ":").c_str(), val.c_str());
}

void PrintTime(const std::string &s, double seconds) {
  if (shmem_my_pe() == 0) 
    printf("PE: %d | %-21s%3.5lf\n", shmem_my_pe(), (s + ":").c_str(), seconds);
}

void PrintStep(const std::string &s, int64_t count) {
  if (shmem_my_pe() == 0) 
    printf("PE: %d | %-14s%14" PRId64 "\n", shmem_my_pe(), (s + ":").c_str(), count);
}

void PrintStep(int step, double seconds, int64_t count = -1) {
  if (shmem_my_pe() == 0) {
    if (count != -1)
      printf("PE: %d | %5d%11" PRId64 "  %10.5lf\n", shmem_my_pe(), step, count, seconds);
    else
      printf("PE: %d | %5d%23.5lf\n", shmem_my_pe(), step, seconds);
  }
}

void PrintStep(const std::string &s, double seconds, int64_t count = -1) {
  if (shmem_my_pe() == 0) {
    if (count != -1)
      printf("PE: %d | %5s%11" PRId64 "  %10.5lf\n", shmem_my_pe(), s.c_str(), count, seconds);
    else
      printf("PE: %d | %5s%23.5lf\n", shmem_my_pe(), s.c_str(), seconds);
  }
}

// Runs op and prints the time it took to execute labelled by label
#define TIME_PRINT(label, op) {   \
  Timer t_;                       \
  t_.Start();                     \
  (op);                           \
  t_.Stop();                      \
  PrintTime(label, t_.Seconds()); \
}


template <typename T_>
class RangeIter {
  T_ x_;
 public:
  explicit RangeIter(T_ x) : x_(x) {}
  bool operator!=(RangeIter const& other) const { return x_ != other.x_; }
  T_ const& operator*() const { return x_; }
  RangeIter& operator++() {
    ++x_;
    return *this;
  }
};

template <typename T_>
class Range{
  T_ from_;
  T_ to_;
 public:
  explicit Range(T_ to) : from_(0), to_(to) {}
  Range(T_ from, T_ to) : from_(from), to_(to) {}
  RangeIter<T_> begin() const { return RangeIter<T_>(from_); }
  RangeIter<T_> end() const { return RangeIter<T_>(to_); }
};

// Return the number of elements that would be allocated to PE
// p with NPES n using round robin partitioning
template <typename T_>
T_ round_robin_size(T_ num_edges, int p, int n) {
  T_ k = num_edges - (num_edges % n);
  T_ I = k / n;
  T_ r = num_edges - k;
  T_ counts[n];
  for (int i = 0; i < n; i++)
    counts[i] = I;
  int j = 0;
  while (j < r) {
    counts[j]++;
    j++;
  }
  return counts[p];  
}


#endif  // UTIL_H_
