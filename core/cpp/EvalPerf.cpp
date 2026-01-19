#include <chrono>      // Standard library for high-resolution timing (steady_clock, duration)
#include <iostream>    // Input/Output stream for printing results to console
#include <x86intrin.h> // Provides access to x86 intrinsic instructions (e.g., RDTSC for CPU cycles)

#include "EvalPerf.hpp"

/**
 * Reads the Time Stamp Counter (TSC) of the CPU.
 * * This function uses inline assembly to execute the 'rdtsc' instruction,
 * which returns the number of clock cycles since the last processor reset.
 * * @return uint64_t The 64-bit unsigned integer representing CPU cycles.
 */
uint64_t rdtsc()
{
      unsigned int lo, hi;

      // Execute rdtsc instruction:
      // "=a" (lo) : stores the lower 32 bits of the TSC in the 'eax' register
      // "=d" (hi) : stores the upper 32 bits of the TSC in the 'edx' register
      __asm__ __volatile__("rdtsc" : "=a"(lo), "=d"(hi));

      // Combine the two 32-bit values into a single 64-bit result
      return ((uint64_t)hi << 32) | lo;
}

// Constructor: Default initialization
EvalPerf::EvalPerf() {}

/**
 * Starts the performance monitoring.
 * Captures the current wall-clock time and CPU cycle count simultaneously.
 */
void EvalPerf::start()
{
      _t0 = std::chrono::high_resolution_clock::now(); // Capture system time
      _c0 = rdtsc();                                   // Capture hardware cycles
}

/**
 * Stops the performance monitoring.
 * Captures the ending wall-clock time and CPU cycle count.
 */
void EvalPerf::stop()
{
      _t1 = std::chrono::high_resolution_clock::now(); // Capture system time
      _c1 = rdtsc();                                   // Capture hardware cycles
}

/**
 * Resets the measurement by aligning start and end points.
 * This effectively sets the duration and cycle delta to zero.
 */
void EvalPerf::clear()
{
      _t0 = _t1;
      _c1 = _c0;
}

/**
 * Returns the elapsed time in milliseconds.
 * Uses std::chrono duration conversion to transform the time_point delta.
 */
double EvalPerf::milliseconds() const
{
      return std::chrono::duration<double, std::milli>(_t1 - _t0).count();
}

/**
 * Returns the elapsed time in seconds.
 * Default duration unit for std::chrono is seconds.
 */
double EvalPerf::seconds() const
{
      return std::chrono::duration<double>(_t1 - _t0).count();
}

/**
 * Returns the total number of CPU cycles elapsed between start() and stop().
 */
uint64_t EvalPerf::cycles() const
{
      return _c1 - _c0;
}

/**
 * Calculates Cycles Per Instruction (CPI).
 * Represents the average number of clock cycles required to execute one operation.
 * @param N The total number of operations performed.
 */
double EvalPerf::CPI(size_t N) const
{
      return double(_c1 - _c0) / N;
}

/**
 * Calculates Instructions Per Cycle (IPC).
 * Measures how many operations are executed in a single clock cycle (Efficiency).
 * @param N The total number of operations performed.
 */
double EvalPerf::IPC(size_t N) const
{
      return N / double(_c1 - _c0);
}

/**
 * Calculates Giga Floating-Point Operations Per Second (Gflops).
 * Measures the computational throughput of the code.
 * @param N The total number of floating-point operations.
 */
double EvalPerf::Gflops(size_t N) const
{
      // Divide by 1e9 to convert from Flops to Gflops
      return double(N) / seconds();
}
