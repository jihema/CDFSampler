/*
 * \file cdf_tester.h
 */

#pragma once

#include <stddef.h>
#include <random>

namespace cdf_sampler
{

namespace test
{

template<typename Sampler> class CDFTester
{
public:
    using Scalar = typename Sampler::Scalar;

    CDFTester(std::random_device::result_type const seed) :
            m_seed(seed), m_sampler(get_random_sampler())
    {
    }

    static bool test_random();

private:

    std::random_device::result_type m_seed;
    Sampler const m_sampler;

    Sampler get_random_sampler();

    static std::string get_sampler_name();

    /**
     * \brief Uses sampler to compute the MC integral of a function identically equal to 1.
     *
     * The result should be approximately equal to the reference value which is the area of the integration domain.
     */
    void check_sum_one(Sampler const& sampler, typename Sampler::Domain const& domain, size_t const num_samples);


    void check_sum_sin(Sampler const& sampler, typename Sampler::Domain const& domain, size_t const num_samples);


    static std::vector<Scalar> wiggle(Scalar const x_min, Scalar const x_max, size_t const x_size, std::mt19937& random)
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
