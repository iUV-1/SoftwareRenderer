// stolen from a CS147 project

#pragma once

#include <chrono>
#include <cstdint>

class CycleTimer {
public:
    using SysClock = std::uint64_t;
    using Clock = std::chrono::steady_clock;
    using Seconds = std::chrono::duration<double>;

    // Returns seconds elapsed since the first call, using a monotonic clock.
    // Thread-safe since C++11 (static local initialization is thread-safe).
    static inline double currentSeconds() {
        static const Clock::time_point t0 = Clock::now(); // ran exactly once, at the start of the program since its static
        return Seconds{Clock::now() - t0}.count();
    }

    // Utility class: no instances.
    CycleTimer() = delete;
    ~CycleTimer() = delete;
    CycleTimer(const CycleTimer&) = delete;
    CycleTimer& operator=(const CycleTimer&) = delete;
};