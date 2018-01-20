/**
 * \file CDFSamplerTest.h
 */

#pragma once

#include <stddef.h>
#include <random>
#include <vector>


namespace cdf_sampler
{

namespace test
{

template<typename Scalar>
class CDFSamplerTest
{
public:
    static void test_cdf_sampler_1D();
    static void test_cdf_sampler_1D_again();
    static void test_cdf_sampler_2D();

private:
    template<typename Sampler>
    static void check_sum(Sampler const& sampler,
            typename Sampler::Domain const& domain, size_t const num_samples,
            unsigned int const* const seed = nullptr);

    static std::vector<Scalar> wiggle(Scalar const x_min, Scalar const x_max,
            size_t const x_size, std::mt19937& random)
    {
        std::vector<Scalar> x(x_size);
        std::uniform_real_distribution<Scalar> uniform_01;
        for (size_t i = 0; i < x_size; ++i)
        {
            Scalar const step = (x_max - x_min) / (x_size - 1);
            x[i] = x_min + i * step;
            if (i > 0 && i < x_size - 1)
            {
                x[i] += (uniform_01(random) - 0.5) * step;
            }
        }
        return x;
    }
};

}

}
