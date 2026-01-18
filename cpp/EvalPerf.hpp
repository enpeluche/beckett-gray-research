#ifndef EVALPERF_HPP
#define EVALPERF_HPP
#include <chrono>
#include <iostream>
#include <x86intrin.h>

class EvalPerf
{
    public:
        EvalPerf();
        void start();
        void stop();
        void clear();
        double milliseconds() const;
        double seconds() const;
        uint64_t cycles() const;
        double CPI(size_t N) const;
        double IPC(size_t N) const;
        double Gflops(size_t N) const;

    protected:

    private:
        std::chrono::time_point<std::chrono::high_resolution_clock> _t0, _t1;
        uint64_t _c0, _c1;
};

#endif // EVALPERF_HPP
