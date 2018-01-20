/**
 * \file RectifiedSamplingTest.h
 */

#pragma once

#include <iostream>
#include <cmath>
#include "../src/vec.h"

namespace cdf_sampler
{

namespace test
{

template<typename Scalar>
class RectifiedSamplingTest
{
public:
    static void test_cdf_sampler_1D();
    static void test_cdf_sampler_1D_again();
    static void test_cdf_sampler_2D();

private:
    template<typename Sampler>
    static void check_sum(Sampler const& sampler,
            typename Sampler::Domain const& domain, size_t const num_samples =
                    10000);
};

}

}
