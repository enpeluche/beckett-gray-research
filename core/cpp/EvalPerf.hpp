#ifndef EVALPERF_HPP
#define EVALPERF_HPP

#include <chrono>      // Standard library for high-resolution timing (steady_clock, duration)
#include <iostream>    // Input/Output stream for printing results to console
#include <x86intrin.h> // Provides access to x86 intrinsic instructions (e.g., RDTSC for CPU cycles)

/**
 * Class to measure code performance using hardware counters and high-resolution timers.
 */
class EvalPerf
{
public:
    EvalPerf(); // Constructor: initializes the performance monitor

    // --- Measurement Control ---
    void start(); // Record current time and CPU cycle count
    void stop();  // Record end time and CPU cycle count
    void clear(); // Reset all internal counters

    // --- Time Metrics ---
    double milliseconds() const; // Elapsed time in milliseconds
    double seconds() const;      // Elapsed time in seconds

    // --- CPU Performance Metrics ---
    uint64_t cycles() const; // Total CPU cycles consumed between start/stop

    /**
     * CPI (Cycles Per Instruction): Average cycles to execute one instruction.
     * @param N Number of instructions (or operations) processed.
     */
    double CPI(size_t N) const;

    /**
     * IPC (Instructions Per Cycle): Reciprocal of CPI.
     * Measures instruction-level parallelism.
     */
    double IPC(size_t N) const;

    /**
     * Gflops: Billion floating-point operations per second.
     * Measures computational throughput.
     */
    double Gflops(size_t N) const;

protected:
private:
    // High-resolution timestamps for wall-clock time measurement
    std::chrono::time_point<std::chrono::high_resolution_clock> _t0, _t1;

    // CPU cycle counters (usually retrieved via RDTSC instruction)
    uint64_t _c0, _c1;
};

#endif // EVALPERF_HPP
