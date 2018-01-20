/**
 * \file monotone_cubic_interpolator.h
 */

#pragma once

#include <vector>

namespace dneg
{

/**
 * \brief Piecewise monotone cubic Hermite interpolation.
 *
 * Implements the Fritsch-Carlson algorithm (https://en.wikipedia.org/wiki/Monotone_cubic_interpolation).
 * Interpolation is guaranteed to be \f$ C^1 \f$ and monotone on any interval where the input data is monotone.
 */
template<typename Scalar>
class MonotoneCubicInterpolator
{
public:

    MonotoneCubicInterpolator(const std::vector<Scalar> & x, const std::vector<Scalar> & f);

    Scalar evaluate(Scalar x) const;

    inline std::vector<Scalar>& get_fp()
    {
        return m_fp;
    }

private:
    std::vector<Scalar> m_x; ///< Abscissa.
    std::vector<Scalar> m_f; ///< Values.
    std::vector<Scalar> m_fp; ///< Derivatives.

    inline static bool different_sign(Scalar a, Scalar b)
    {
        return (a < 0 && b > 0) || (a > 0 && b < 0);
    }

};

}
