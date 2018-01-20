/*
 * \file vec.h
 */

#pragma once

#include <array>
#include <ostream>

namespace cdf_sampler
{

template<typename Scalar> using vec2 = std::array<Scalar, 2>;
template<typename Scalar> using vec3 = std::array<Scalar, 3>;
template<typename Scalar> using box2 = std::array<vec2<Scalar>, 2>;

using vec2f = vec2<float>;
using vec2d = vec2<double>;
using box2f = box2<float>;
using box2d = box2<double>;

template<typename Scalar>
std::ostream& operator<<(std::ostream& os, vec2<Scalar> const& x)
{
    return os << "{" << x[0] << ", " << x[1] << "}";
}

template<typename Scalar>
vec2<Scalar> operator+(vec2<Scalar> const& x1, vec2<Scalar> const& x2)
{
    return vec2<Scalar> { x1[0] + x2[0], x1[1] + x2[2] };
}

template<typename Scalar0, typename Scalar>
vec2<Scalar> operator*(Scalar0 const s, vec2<Scalar> const& x)
{
    return vec2<Scalar> { (Scalar) s * x[0], (Scalar) s * x[2] };
}

template<typename Domain, typename Scalar>
Scalar compute_area(Domain const& domain);

template<>
inline float compute_area(vec2f const& domain)
{
    return domain[1] - domain[0];
}

template<>
inline double compute_area(vec2d const& domain)
{
    return domain[1] - domain[0];
}

template<>
inline float compute_area(box2f const& domain)
{
    return (domain[1][0] - domain[0][0]) * (domain[1][1] - domain[0][1]);
}

template<>
inline double compute_area(box2d const& domain)
{
    return (domain[1][0] - domain[0][0]) * (domain[1][1] - domain[0][1]);
}

template<typename Domain>
Domain intersect(Domain const& domain1, Domain const& domain2);

template<typename Scalar>
inline vec2<Scalar> intersect(vec2<Scalar> const& domain1,
        vec2<Scalar> const& domain2)
{
    return vec2<Scalar> { std::max(domain1[0], domain2[0]), std::min(domain1[1],
            domain2[1]) };
}

template<typename Scalar>
inline box2<Scalar> intersect(box2<Scalar> const& domain1,
        box2<Scalar> const domain2)
{
    return box2<Scalar> { std::max(domain1[0][0], domain2[0][0]), std::max(
            domain1[0][1], domain2[0][1]), std::min(domain1[1][0],
            domain2[1][0]), std::min(domain1[1][1], domain2[1][1]) };
}

}
