/*
 * \file scoped_timer.h
 */

#pragma once

#include <chrono>
#include <functional>

namespace cdf_sampler
{

namespace test
{

class ScopedTimer
{

public:

    ScopedTimer(std::function<void(int64_t)> callback) :
            m_start_time(std::chrono::high_resolution_clock::now()), //
            m_callback(callback)
    {
    }

    ~ScopedTimer()
    {
        std::chrono::high_resolution_clock::time_point const now_time =
                std::chrono::high_resolution_clock::now();
        int64_t const duration = std::chrono::duration_cast<
                std::chrono::nanoseconds>(now_time - m_start_time).count();
        m_callback(duration);
    }

private:

    std::chrono::high_resolution_clock::time_point const m_start_time;
    std::function<void(int64_t)> const m_callback;
};

}

}
