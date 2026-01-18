#include <chrono>
#include <iostream>
#include <x86intrin.h>
#include "EvalPerf.hpp"
uint64_t rdtsc(){
    unsigned int lo,hi;
    __asm__ __volatile__ ("rdtsc" : "=a" (lo), "=d" (hi));
    return ((uint64_t)hi << 32) | lo;
}

EvalPerf::EvalPerf() {};

void EvalPerf::start(){
      _t0 = std::chrono::high_resolution_clock::now();
      _c0=rdtsc();
      }


void EvalPerf::stop() {
      _t1 = std::chrono::high_resolution_clock::now();
      _c1=rdtsc();
      }


void EvalPerf::clear() {
      _t0=_t1;
      _c1=_c0;
      }

double EvalPerf::milliseconds() const {
      return std::chrono::duration<double, std::milli>(_t1-_t0).count();
      }

double EvalPerf::seconds() const {
      return std::chrono::duration<double>(_t1-_t0).count();
      }

uint64_t EvalPerf::cycles() const {
      return _c1-_c0;
      }

double EvalPerf::CPI(size_t N) const {
      return double(_c1-_c0)/N;
      }

double EvalPerf::IPC(size_t N) const {
      return N/double(_c1-_c0);
      }

double EvalPerf::Gflops(size_t N) const {
      return double(N)/seconds();
      }



