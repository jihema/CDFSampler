/**
 * \file RectifiedSamplingTest.h
 */

#pragma once

#include <iostream>
#include <cmath>

#include <gmath_vec2.h>
#include <gmath_bbox2.h>

namespace dneg
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
    static void test_ies_interpolation();

private:
    template<typename Sampler>
    static void check_sum(Sampler const& cdf_sampler, typename Sampler::Domain const& domain);
};

template<typename Domain, typename Scalar>
Scalar compute_area(Domain const& domain);

template<>
inline float compute_area(GMathVec2f const& domain)
{
    return domain[1] - domain[0];
}

template<>
inline double compute_area(GMathVec2d const& domain)
{
    return domain[1] - domain[0];
}

template<>
inline float compute_area(GMathBbox2f const& domain)
{
    return domain.compute_area();
}

template<>
inline double compute_area(GMathBbox2d const& domain)
{
    return domain.compute_area();
}

template<typename Domain>
Domain intersect(Domain const& domain1, Domain const& domain2);

template<typename Scalar>
inline GMathVec2<Scalar> intersect(GMathVec2<Scalar> const& domain1,
        GMathVec2<Scalar> const& domain2)
{
    return GMathVec2<Scalar>(std::max(domain1[0], domain2[0]), std::min(domain1[1], domain2[1]));
}

template<typename Scalar>
inline GMathBbox2<Scalar> intersect(GMathBbox2<Scalar> const& domain1,
        GMathBbox2<Scalar> const domain2)
{
    return const_cast<GMathBbox2<Scalar>&>(domain1).intersects(domain2);
}

}

} /* namespace dneg */
