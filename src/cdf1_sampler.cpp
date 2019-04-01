/**
 * \file cdf1_sampler.cpp
 */

#include "cdf1_sampler.h"

#include <stddef.h>

namespace cdf_sampler
{

template<typename Scalar>
vec2<Scalar> const CDF1Sampler<Scalar>::max_domain {
        -std::numeric_limits<Scalar>::max(), std::numeric_limits<Scalar>::max() };

template<typename Scalar>
void CDF1SamplerUniform<Scalar>::init(Domain const& x_range,
        std::vector<Scalar> const& f)
{
    m_x_min = x_range[0];
    m_x_max = x_range[1];

    this->m_pdf = f; // To be normalized.
    this->m_cdf.resize(f.size());

    Scalar const step = (m_x_max - m_x_min) / (this->m_pdf.size() - 1);

    Scalar previous_cdf = 0.;
    for (size_t xi = 0; xi < this->m_pdf.size(); ++xi)
    {
        Scalar delta_x;
        if (xi == 0 || xi == this->m_pdf.size() - 1)
        {
            delta_x = 0.5 * step;
        }
        else
        {
            delta_x = step;
        }

        previous_cdf = CDF1Sampler<Scalar>::get_cdf(xi) = previous_cdf
                + delta_x * this->get_pdf(xi);
    }

    this->normalize();
}

template<typename Scalar>
Scalar CDF1SamplerUniform<Scalar>::sample(Scalar const in_s, Scalar& pdf,
        Domain const& domain) const
{
    assert(0. <= in_s && in_s <= 1.);

    Scalar const domain_min = std::max(m_x_min, domain[0]);
    size_t ix_min = -1;
    Scalar theta_x_min = std::numeric_limits<Scalar>::signaling_NaN();
    Scalar const cdf0 = get_cdf(domain_min, ix_min, theta_x_min);

    Scalar const domain_max = std::min(m_x_max, domain[1]);
    size_t ix_max = -1;
    Scalar theta_x_max = std::numeric_limits<Scalar>::signaling_NaN();
    Scalar const cdf1 = get_cdf(domain_max, ix_max, theta_x_max);

    if (cdf1 - cdf0 < std::numeric_limits<Scalar>::epsilon())
    {
        pdf = 0.;
        return 0.5 * (domain_min + domain_max);
    }

    if (ix_min == ix_max) // First case: both x_min and x_max are in the same interval.
                          // This means we are sampling uniformly between domain_min and domain_max.
    {
        pdf = 1. / (domain_max - domain_min);
        return linear_interpolation(domain_min, domain_max, in_s);
    }
    else // Second case: [x_min, x_max] straddle at least one point of the interleaved x range.
    {
        // Let F(x) := \int_{x_min}^{x_max} pdf(t, y) dt.
        // We create a temporary cdf values range, starting with 0 = F(x_min), ending with F(x_max),
        // and with all the interleaved x range values in between.
        std::vector<Scalar> cdfx;
        cdfx.reserve(ix_max - ix_min + 2);
        cdfx.push_back(0.);
        for (auto xi = ix_min; xi < ix_max; ++xi)
        {
            Scalar const cdfxy = this->m_cdf[xi];
            cdfx.push_back(cdfxy - cdf0);
        }
        cdfx.push_back(cdf1 - cdf0);
        Interpol<std::vector<Scalar>> const xax = find_range(in_s * cdfx.back(),
                cdfx);

        // Convert the interpol to one on the cdf positions vector.
        size_t hit_ix;
        Scalar hit_theta;
        if (xax.m_first == cdfx.begin())
        {
            hit_ix = ix_min;
            hit_theta = theta_x_min + xax.m_theta * (1. - theta_x_min);
        }
        else if (xax.m_first == cdfx.end() - 2)
        {
            hit_ix = ix_max;
            hit_theta = xax.m_theta * theta_x_max;
        }
        else
        {
            hit_ix = ix_min + (xax.m_first - cdfx.begin());
            hit_theta = xax.m_theta;
        }

        pdf = this->m_pdf[hit_ix] / (cdf1 - cdf0);

        return linear_interpolation(interleaved_x(hit_ix),
                interleaved_x(hit_ix + 1), hit_theta);
    }
}

template<typename Scalar>
void CDF1SamplerNonUniform<Scalar>::init(std::vector<Scalar> const& x,
        std::vector<Scalar> const& f)
{
    assert(f.size() == x.size());

    m_x = x;
    this->m_pdf = f; // To be normalized.
    this->m_cdf.resize(f.size());

    m_interleaved_x = CDFSampler<Scalar>::interleave(m_x);

    Scalar previous_cdf = 0.;

    for (size_t xi = 0; xi < m_x.size(); ++xi)
    {
        Scalar delta_x;
        if (xi == 0)
        {
            delta_x = 0.5 * (m_x[xi + 1] - m_x[xi]);
        }
        else if (xi == m_x.size() - 1)
        {
            delta_x = 0.5 * (m_x[xi] - m_x[xi - 1]);
        }
        else
        {
            delta_x = 0.5 * (m_x[xi + 1] - m_x[xi - 1]);
        }

        previous_cdf = CDF1Sampler<Scalar>::get_cdf(xi) = previous_cdf
                + delta_x * CDF1Sampler<Scalar>::get_pdf(xi);
    }

    this->normalize();
}

template<typename Scalar>
Scalar CDF1SamplerNonUniform<Scalar>::sample(Scalar const in_s, Scalar& pdf,
        vec2<Scalar> const& domain) const
{
    assert(0. <= in_s && in_s <= 1.);

    // TODO: we could cache the following, depending only on the domain.
    Interpol<std::vector<Scalar>> const x_min = find_range(domain[0],
            m_interleaved_x);
    Interpol<std::vector<Scalar>> const x_max = find_range(domain[1],
            m_interleaved_x);

    Scalar const cdf0 = get_cdf(x_min);
    Scalar const cdf1 = get_cdf(x_max);

    if (cdf1 - cdf0 < std::numeric_limits<Scalar>::epsilon())
    {
        pdf = 0.;
        return 0.5 * (domain[0] + domain[1]);
    }

    Interpol<std::vector<Scalar>> hit_x; // This will receive an interpol for the sample in the interleaved x range.

    if (x_min.m_first == x_max.m_first) // First case: both x_min and x_max are in the same interval.
    {
        hit_x.m_first = x_min.m_first;
        hit_x.m_theta = x_min.m_theta + in_s * (x_max.m_theta - x_min.m_theta);
    }
    else // Second case: [x_min, x_max] straddle at least one point of the interleaved x range.
    {
        // Let F(x) := \int_{x_min}^{x_max} pdf(t, y) dt.
        // We create a temporary cdf values range, starting with 0 = F(x_min), ending with F(x_max),
        // and with all the interleaved x range values in between.
        std::vector<Scalar> cdfx;
        cdfx.reserve(x_max.m_first - x_min.m_first + 2);
        cdfx.push_back(0.);
        for (auto xi = x_min.m_first + 1; xi <= x_max.m_first; ++xi)
        {
            Scalar const cdfxy = this->m_cdf[xi - m_interleaved_x.begin() - 1];
            cdfx.push_back(cdfxy - cdf0);
        }
        cdfx.push_back(cdf1 - cdf0);
        Interpol<std::vector<Scalar>> const xax = find_range(in_s * cdfx.back(),
                cdfx);

        // Convert the interpol to one on the cdf positions vector.
        if (xax.m_first == cdfx.begin()) // Hit below first x-range position.
        {
            hit_x.m_first = x_min.m_first;
            hit_x.m_theta = x_min.m_theta + xax.m_theta * (1 - x_min.m_theta);
        }
        else if (xax.m_first == cdfx.end() - 2)
        {
            hit_x.m_first = x_max.m_first;
            hit_x.m_theta = xax.m_theta * x_max.m_theta;
        }
        else
        {
            hit_x.m_first = x_min.m_first + (xax.m_first - cdfx.begin());
            hit_x.m_theta = xax.m_theta;
        }
    }

    pdf = get_pdf(hit_x) / (cdf1 - cdf0);

    return hit_x.value();
}

template class CDF1Sampler<float> ;
template class CDF1Sampler<double> ;
template class CDF1SamplerUniform<float> ;
template class CDF1SamplerUniform<double> ;
template class CDF1SamplerNonUniform<float> ;
template class CDF1SamplerNonUniform<double> ;

}
