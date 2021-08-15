// Copyright (c) 2015, The Regents of the University of California (Regents)
// See LICENSE.txt for license details

#ifndef TIMER_H_
#define TIMER_H_

#include <chrono>


/*
GAP Benchmark Suite
Class:  Timer
Authors: Scott Beamer, Michael Sutton

Simple timer that wraps std::chrono
*/

/*
DMM-GAPBS
Author: Zach Hansen
Adaptation Notes:
 - All PEs must participate in timer starting/stopping
 - Typically PE 0 reports the time
*/

class Timer {
 public:
  Timer() {}

  void Start() {
    shmem_barrier_all();                                // all pes must be present to start timer
    elapsed_time_ = start_time_ = std::chrono::high_resolution_clock::now();    // what if exec times for this line vary? will each pe have a different start time?
    shmem_barrier_all();
  }

  void Stop() {
    shmem_barrier_all();                                // all pes must arrive to stop timer
    elapsed_time_ = std::chrono::high_resolution_clock::now();
    shmem_barrier_all();
  }

  double Seconds() const {
    return std::chrono::duration_cast<std::chrono::duration<double>>(elapsed_time_ - start_time_).count();
  }

  double Millisecs() const {
    return std::chrono::duration_cast<std::chrono::duration<double, std::milli>>(elapsed_time_ - start_time_).count();
  }

  double Microsecs() const {
    return std::chrono::duration_cast<std::chrono::duration<double, std::micro>>(elapsed_time_ - start_time_).count();
  }

 private:
  std::chrono::high_resolution_clock::time_point start_time_, elapsed_time_;
};

// Times op's execution using the timer t
#define TIME_OP(t, op) { t.Start(); (op); t.Stop(); }

#endif  // TIMER_H_
