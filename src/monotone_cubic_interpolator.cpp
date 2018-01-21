/**
 * \file monotone_cubic_interpolator.cpp
 */

#include "monotone_cubic_interpolator.h"
#include "interpol.h"
#include <limits>
#include <cmath>
#include <assert.h>

using namespace std;

namespace cdf_sampler
{

template<typename Scalar>
MonotoneCubicInterpolator<Scalar>::MonotoneCubicInterpolator(
        const vector<Scalar> & x, const vector<Scalar> & f) :
        m_x(x), m_f(f)
{
    assert(m_x.size() == m_f.size());
    assert(m_x.size() >= 2);

    // First estimate of derivatives. First and last points get right and left finite difference
    // respectively; intermediary points get average of both sides.
    m_fp.reserve(m_f.size());
    typename std::vector<Scalar>::const_iterator x_it = m_x.begin();
    typename std::vector<Scalar>::const_iterator f_it = m_f.begin();
    Scalar delta = (*f_it - *(f_it + 1)) / (*x_it - *(x_it + 1));
    m_fp.push_back(delta);
    for (; x_it != m_x.end() - 2;)
    {
        m_fp.push_back(delta);
        x_it++;
        f_it++;
        delta = (*f_it - *(f_it + 1)) / (*x_it - *(x_it + 1));
        m_fp.back() += delta;
        m_fp.back() *= 0.5;
    }
    m_fp.push_back(delta);

    // Adjust derivatives for monotonicity (Fritsch-Carlson algorithm).
    for (size_t i = 0; i < m_x.size() - 1; ++i)
    {
        Scalar const delta = (m_f[i + 1] - m_f[i]) / (m_x[i + 1] - m_x[i]);
        if (fabs(delta) < std::numeric_limits<Scalar>::epsilon()) // Locally constant between this point and the next.
        {
            m_fp[i] = m_fp[i + 1] = 0.;
        }
        else if (different_sign(m_fp[i], delta)) // This point is a local extremum.
        {
            m_fp[i] = 0.;
        }
        else if (different_sign(m_fp[i + 1], delta)) // Next one is a local extremum.
        {
            m_fp[i + 1] = 0.;
        }
        else
        {
            Scalar const alpha = m_fp[i] / delta;
            Scalar const beta = m_fp[i + 1] / delta;
            Scalar const sqnorm = alpha * alpha + beta * beta;
            if (sqnorm > 9.)
            {
                Scalar const tau = 3. / sqrt(sqnorm);
                m_fp[i] = tau * alpha * delta;
                m_fp[i + 1] = tau * beta * delta;
            }
        }
    }
}

template<typename Scalar>
Scalar MonotoneCubicInterpolator<Scalar>::evaluate(Scalar x) const
{
    Interpol<std::vector<Scalar>> const x_it = find_range(x, m_x);
    size_t const x_idx = x_it.m_first - m_x.begin();
    return cubic_interpolation(m_f[x_idx], m_f[x_idx + 1], m_fp[x_idx],
            m_fp[x_idx + 1], x_it.m_theta);
}

template class MonotoneCubicInterpolator<float> ;
template class MonotoneCubicInterpolator<double> ;

}
