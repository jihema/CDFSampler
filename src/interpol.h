/**
 * \file interpol.h
 */

#pragma once

#include <assert.h>
#include <ostream>

namespace cdf_sampler
{

/**
 * \brief Represents an interpolated position in a container.
 *
 * The value at that position can be computed as
 * \code{*m_first + (*(m_first + 1)- * m_first) * m_theta}\endcode.
 */
template<typename Container>
class Interpol
{
public:
    Interpol()
    {
    }

    Interpol(typename Container::const_iterator const& first, typename Container::value_type theta) :
            m_first(first), m_theta(theta)
    {
    }

    typename Container::value_type value() const
    {
        return *m_first + (*(m_first + 1) - *m_first) * m_theta;
    }

    typename Container::const_iterator m_first;
    typename Container::value_type m_theta;
};

/**
 * \brief Finds interpolating range in an ordered container.
 *
 * The container must have bidirectional iterator.
 *
 * \return A pair (iterator, theta) with 0 <= theta <= 1 and such that
 * x = (1 - theta) (*iterator) + theta (*(iterator+1)).
 * In particular, iterator >= values.begin() and iterator < values.end()-1;
 */
template<typename Container>
static Interpol<Container> find_range(typename Container::value_type const& x,
        Container const& values)
{
    assert(values.size() >= 2);

    using Iterator = typename Container::const_iterator;

    // Now find the brackets for interpolation.
    Iterator const top = std::lower_bound(values.begin(), values.end(), x);
    if(top == values.end()) // Out of range, we'll clamp down to the last value.
    {
        return Interpol<Container>(values.end() - 2, 1.);
    }
    else if(top == values.begin()) // At or below beginning, clamp up to first value.
    {
        return Interpol<Container>(values.begin(), 0.);
    }
    else // In bona fide range between top - 1 and top.
    {
        Iterator const bottom = top - 1;
        return Interpol<Container>(bottom, (x - *bottom) / (*top - *bottom));
    }
}

template<typename Scalar, typename Interpolated>
inline Scalar linear_interpolation(Interpolated const& z0, Interpolated const& z1, Scalar theta)
{
    return z0 + (z1 - z0) * theta;
}

template<typename Scalar, typename Interpolated>
inline Scalar cubic_interpolation(Interpolated const& z0, Interpolated const& z1,
        Interpolated const& zp0, Interpolated const& zp1, Scalar t)
{
    Scalar const t2 = t * t;
    Scalar const t3 = t * t2;

    Scalar const h00 = 2 * t3 - 3 * t2 + 1.;
    Scalar const h01 = -2 * t3 + 3 * t2;
    Scalar const h10 = t3 - 2 * t2 + t;
    Scalar const h11 = t3 - t2;

    return h00 * z0 + h01 * z1 + h10 * zp0 + h11 * zp1;
}

template<typename Scalar, typename Interpolated>
inline Scalar bilinear_interpolation(Interpolated const& z00, Interpolated const& z10,
        Interpolated const& z01, Interpolated const& z11, Scalar theta_h, Scalar theta_v)
{
    return (z00 + (z10 - z00) * theta_h) * (1 - theta_v) + (z01 + (z11 - z01) * theta_h) * theta_v;
}

}

