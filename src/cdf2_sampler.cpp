/**
 * \file cdf2_sampler.cpp
 */

#include "cdf2_sampler.h"

#include <stddef.h>
#include <cassert>
#include <limits>

namespace cdf_sampler
{

template<typename Scalar>
void CDF2SamplerUniform<Scalar>::init(Domain const& xy_range, size_t const x_size,
        size_t const y_size, std::vector<Scalar> const& f)
{
    assert(f.size() == x_size * y_size);

    m_x_min = xy_range[0][0];
    m_y_min = xy_range[0][1];
    m_x_max = xy_range[1][0];
    m_y_max = xy_range[1][1];

    this->m_x_size = x_size;
    this->m_y_size = y_size;
    this->m_pdf = f; // To be normalized.
    this->m_cdf.resize(f.size());
    this->m_sum = 0.;

    Scalar const step_x = (m_x_max - m_x_min) / (x_size - 1);
    Scalar const step_y = (m_y_max - m_y_min) / (y_size - 1);

    // First iy = 0.
    {
        Scalar previous_cdf = 0.;
        Scalar const delta_y = 0.5 * step_y;
        for (size_t ix = 0; ix < x_size; ++ix)
        {
            Scalar delta_x;
            if(ix == 0 || ix == x_size - 1)
            {
                delta_x = 0.5 * step_x;
            }
            else
            {
                delta_x = step_x;
            }

            previous_cdf = static_cast<CDF2Sampler<Scalar>*>(this)->get_cdf(ix, 0) = previous_cdf
                    + delta_x * delta_y
                            * static_cast<CDF2Sampler<Scalar> const* const >(this)->get_pdf(ix, 0);
        }
    }

    for (size_t iy = 1; iy < y_size; ++iy)
    {
        Scalar delta_y;
        if(iy == y_size - 1)
        {
            delta_y = 0.5 * step_y;
        }
        else
        {
            delta_y = step_y;
        }

        Scalar previous_cdf_x = 0;
        for (size_t ix = 0; ix < x_size; ++ix)
        {
            Scalar delta_x;
            if(ix == 0 || ix == x_size - 1)
            {
                delta_x = 0.5 * step_x;
            }
            else
            {
                delta_x = step_x;
            }

            previous_cdf_x += delta_x * delta_y
                    * static_cast<CDF2Sampler<Scalar>*>(this)->get_pdf(ix, iy);
            static_cast<CDF2Sampler<Scalar>*>(this)->get_cdf(ix, iy) = static_cast<CDF2Sampler<
                    Scalar>*>(this)->get_cdf(ix, iy - 1) + previous_cdf_x;
        }
    }

    this->normalize();
}

template<typename Scalar>
vec2<Scalar> CDF2SamplerUniform<Scalar>::sample(vec2<Scalar> const& sxy, Scalar& pdf,
        box2<Scalar> const& domain) const
{
    Scalar const sx = sxy[0];
    Scalar const sy = sxy[1];

    assert(0. <= sx && sx <= 1.);
    assert(0. <= sy && sy <= 1.);

    vec2<Scalar> result;

    Scalar const domain_x_min = std::max(m_x_min, domain[0][0]);
    size_t ix_min = -1;
    Scalar theta_x_min = std::numeric_limits<Scalar>::signaling_NaN();

    Scalar const domain_x_max = std::min(m_x_max, domain[1][0]);
    size_t ix_max = -1;
    Scalar theta_x_max = std::numeric_limits<Scalar>::signaling_NaN();

    Scalar const domain_y_min = std::max(m_y_min, domain[0][1]);
    size_t iy_min = -1;
    Scalar theta_y_min = std::numeric_limits<Scalar>::signaling_NaN();

    Scalar const domain_y_max = std::min(m_y_max, domain[1][1]);
    size_t iy_max = -1;
    Scalar theta_y_max = std::numeric_limits<Scalar>::signaling_NaN();

    Scalar const cdf00 = get_cdf(domain_x_min, domain_y_min, ix_min, theta_x_min, iy_min,
            theta_y_min);
    Scalar const cdf10 = get_cdf(domain_x_max, domain_y_min, ix_max, theta_x_max, iy_min,
            theta_y_min);
    Scalar const cdf01 = get_cdf(domain_x_min, domain_y_max, ix_min, theta_x_min, iy_max,
            theta_y_max);
    Scalar const cdf11 = get_cdf(domain_x_max, domain_y_max, ix_max, theta_x_max, iy_max,
            theta_y_max);

    Scalar const domain_proba = (cdf11 - cdf10) + (cdf00 - cdf01);
    if(domain_proba < std::numeric_limits<Scalar>::epsilon())
    {
        pdf = 0;
        return 0.5 * (domain[0] + domain[1]);
    }

    size_t hit_y_first;
    Scalar hit_y_theta;

    if(iy_min == iy_max) // First case: both y_min and y_max are in the same interval.
    {
        hit_y_first = iy_min;
        hit_y_theta = theta_y_min + sy * (theta_y_max - theta_y_min);
    }
    else // Second case: [y_min, y_max] straddle at least one point of the interleaved x range.
    {
        // Let F(y) := \int_{y_min}^{y} \int_{x_min}^{x_max} pdf(x, t) dx dt.
        // We create a temporary range, starting with 0 = F(y_min), ending with F(y_max),
        // and with all values for interleaved y in between.
        std::vector<Scalar> cdfy;
        cdfy.reserve(iy_max - iy_min + 2);
        cdfy.push_back(0.);
        for (auto iy = iy_min; iy < iy_max; ++iy)
        {
            Scalar const cdf0y = this->get_cdf(ix_min, iy);
            Scalar const cdf1y = this->get_cdf(ix_max, iy);
            cdfy.push_back((cdf1y - cdf10) + (cdf00 - cdf0y));
        }
        cdfy.push_back(domain_proba);

        // Sampling on the y-axis. This gives us an interpol on the extended y-positions vector.
        Interpol<std::vector<Scalar>> const yay = find_range(sy * cdfy.back(), cdfy);

        // Convert the interpol to one on the cdf positions vector.
        if(yay.m_first == cdfy.begin()) // Hit below first y-range position.
        {
            hit_y_first = iy_min;
            hit_y_theta = theta_y_min + yay.m_theta * (1 - theta_y_min);
        }
        else if(yay.m_first == cdfy.end() - 2)
        {
            hit_y_first = iy_max;
            hit_y_theta = yay.m_theta * theta_y_max;
        }
        else
        {
            hit_y_first = iy_min + (yay.m_first - cdfy.begin());
            hit_y_theta = yay.m_theta;
        }
    }

    result[1] = linear_interpolation(interleaved_y(hit_y_first), interleaved_y(hit_y_first + 1),
            hit_y_theta);

    // Still need to sample in x.
    size_t hit_x_first;
    Scalar hit_x_theta;

    if(ix_min == ix_max) // First case: both x_min and x_max are in the same interval.
    {
        hit_x_first = ix_min;
        hit_x_theta = theta_x_min + sx * (theta_x_max - theta_x_min);
    }
    else // Second case: [x_min, x_max] straddle at least one point of the interleaved x range.
    {
        // Let F(x) := \int_{x_min}^{x_max} pdf(t, y) dt (y being the sample obtained above).
        // We create a temporary cdf values range, starting with 0 = F(x_min), ending with F(x_max),
        // and with all values for interleaved x in between.
        std::vector<Scalar> cdfx;
        cdfx.reserve(ix_max - ix_min + 2);
        cdfx.push_back(0.);

        Scalar const cdf0y = band_cdf(ix_min, theta_x_min, hit_y_first);
        for (auto ix = ix_min; ix < ix_max; ++ix)
        {
            Scalar const cdfxy = band_cdf(ix, 0., hit_y_first);
            cdfx.push_back(cdfxy - cdf0y);
        }
        Scalar const cdf1y = band_cdf(ix_max, theta_x_max, hit_y_first);
        cdfx.push_back(cdf1y - cdf0y);

        Interpol<std::vector<Scalar>> const xax = find_range(sx * cdfx.back(), cdfx);

        // Convert the interpol to one on the cdf positions vector.
        if(xax.m_first == cdfx.begin()) // Hit below first x-range position.
        {
            hit_x_first = ix_min;
            hit_x_theta = theta_x_min + xax.m_theta * (1 - theta_x_min);
        }
        else if(xax.m_first == cdfx.end() - 2)
        {
            hit_x_first = ix_max;
            hit_x_theta = xax.m_theta * theta_x_max;
        }
        else
        {
            hit_x_first = ix_min + (xax.m_first - cdfx.begin());
            hit_x_theta = xax.m_theta;
        }
    }

    result[0] = linear_interpolation(interleaved_x(hit_x_first), interleaved_x(hit_x_first + 1),
            hit_x_theta);

    pdf = get_pdf(hit_x_first, hit_y_first) / domain_proba;

    return result;
}

template<typename Scalar>
void CDF2SamplerNonUniform<Scalar>::init(std::vector<Scalar> const& x, std::vector<Scalar> const& y,
        std::vector<Scalar> const& f)
{
    assert(f.size() == x.size() * y.size());

    m_x = x;
    this->m_x_size = m_x.size();
    m_y = y;
    this->m_y_size = m_y.size();
    this->m_pdf = f; // To be normalized.
    this->m_cdf.resize(f.size());
    this->m_sum = 0.;

    m_interleaved_x = CDFSampler<Scalar>::interleave(m_x);
    m_interleaved_y = CDFSampler<Scalar>::interleave(m_y);

    // First iy = 0.
    {
        Scalar previous_cdf = 0;
        Scalar const delta_y = m_y[1] - m_y[0];
        for (size_t ix = 0; ix < m_x.size(); ++ix)
        {
            Scalar delta_x;
            if(ix == 0)
            {
                delta_x = m_x[ix + 1] - m_x[ix];
            }
            else if(ix == m_x.size() - 1)
            {
                delta_x = m_x[ix] - m_x[ix - 1];
            }
            else
            {
                delta_x = m_x[ix + 1] - m_x[ix - 1];
            }

            previous_cdf = static_cast<CDF2Sampler<Scalar>*>(this)->get_cdf(ix, 0) = previous_cdf
                    + 0.25 * delta_x * delta_y
                            * static_cast<CDF2Sampler<Scalar> const* const >(this)->get_pdf(ix, 0);
        }
    }

    for (size_t iy = 1; iy < m_y.size(); ++iy)
    {
        Scalar delta_y;
        if(iy == m_y.size() - 1)
        {
            delta_y = m_y[iy] - m_y[iy - 1];
        }
        else
        {
            delta_y = m_y[iy + 1] - m_y[iy - 1];
        }

        Scalar previous_cdf_x = 0;
        for (size_t ix = 0; ix < m_x.size(); ++ix)
        {
            Scalar delta_x;
            if(ix == 0)
            {
                delta_x = m_x[ix + 1] - m_x[ix];
            }
            else if(ix == m_x.size() - 1)
            {
                delta_x = m_x[ix] - m_x[ix - 1];
            }
            else
            {
                delta_x = m_x[ix + 1] - m_x[ix - 1];
            }

            previous_cdf_x += 0.25 * delta_x * delta_y
                    * static_cast<CDF2Sampler<Scalar>*>(this)->get_pdf(ix, iy);
            static_cast<CDF2Sampler<Scalar>*>(this)->get_cdf(ix, iy) = static_cast<CDF2Sampler<
                    Scalar>*>(this)->get_cdf(ix, iy - 1) + previous_cdf_x;
        }
    }

    this->normalize();
}

template<typename Scalar>
vec2<Scalar> CDF2SamplerNonUniform<Scalar>::sample(vec2<Scalar> const& sxy, Scalar& pdf,
        box2<Scalar> const& domain) const
{
    Scalar const sx = sxy[0];
    Scalar const sy = sxy[1];

    assert(0. <= sx && sx <= 1.);
    assert(0. <= sy && sy <= 1.);

    vec2<Scalar> result;

    // TODO: we could cache the following, depending only on the domain.
    Interpol<std::vector<Scalar>> const x_min = find_range(domain[0][0], m_interleaved_x);
    Interpol<std::vector<Scalar>> const x_max = find_range(domain[1][0], m_interleaved_x);
    Interpol<std::vector<Scalar>> const y_min = find_range(domain[0][1], m_interleaved_y);
    Interpol<std::vector<Scalar>> const y_max = find_range(domain[1][1], m_interleaved_y);

    Scalar const cdf00 = get_cdf(x_min, y_min);
    Scalar const cdf10 = get_cdf(x_max, y_min);
    Scalar const cdf01 = get_cdf(x_min, y_max);
    Scalar const cdf11 = get_cdf(x_max, y_max);

    Scalar const domain_proba = (cdf11 - cdf10) + (cdf00 - cdf01);
    if(domain_proba < std::numeric_limits<Scalar>::epsilon())
    {
        pdf = 0;
        return 0.5 * (domain[0] + domain[1]);
    }

    Interpol<std::vector<Scalar>> hit_y; // This will receive an interpol for the sample in the interleaved y range.

    if(y_min.m_first == y_max.m_first) // First case: both y_min and y_max are in the same interval.
    {
        hit_y.m_first = y_min.m_first;
        hit_y.m_theta = y_min.m_theta + sy * (y_max.m_theta - y_min.m_theta);
    }
    else // Second case: [y_min, y_max] straddle at least one point of the interleaved x range.
    {
        // Let F(y) := \int_{y_min}^{y} \int_{x_min}^{x_max} pdf(x, t) dx dt.
        // We create a temporary range, starting with 0 = F(y_min), ending with F(y_max),
        // and with all values for interleaved y in between.
        std::vector<Scalar> cdfy;
        cdfy.reserve(y_max.m_first - y_min.m_first + 2);
        cdfy.push_back(0.);
        for (auto yi = y_min.m_first + 1; yi <= y_max.m_first; ++yi)
        {
            Interpol<std::vector<Scalar>> const yy(yi, 0.);
            Scalar const cdf0y = get_cdf(x_min, yy);
            Scalar const cdf1y = get_cdf(x_max, yy);
            cdfy.push_back((cdf1y - cdf10) + (cdf00 - cdf0y));
        }
        cdfy.push_back(domain_proba);

        // Sampling on the y-axis. This gives us an interpol on the extended y-positions vector.
        Interpol<std::vector<Scalar>> const yay = find_range(sy * cdfy.back(), cdfy);

        // Convert the interpol to one on the cdf positions vector.
        if(yay.m_first == cdfy.begin()) // Hit below first y-range position.
        {
            hit_y.m_first = y_min.m_first;
            hit_y.m_theta = y_min.m_theta + yay.m_theta * (1 - y_min.m_theta);
        }
        else if(yay.m_first == cdfy.end() - 2)
        {
            hit_y.m_first = y_max.m_first;
            hit_y.m_theta = yay.m_theta * y_max.m_theta;
        }
        else
        {
            hit_y.m_first = y_min.m_first + (yay.m_first - cdfy.begin());
            hit_y.m_theta = yay.m_theta;
        }
    }

    result[1] = hit_y.value();

    // Still need to sample in x.
    Interpol<std::vector<Scalar>> hit_x; // This will receive an interpol for the sample in the interleaved x range.

    if(x_min.m_first == x_max.m_first) // First case: both x_min and x_max are in the same interval.
    {
        hit_x.m_first = x_min.m_first;
        hit_x.m_theta = x_min.m_theta + sx * (x_max.m_theta - x_min.m_theta);
    }
    else // Second case: [x_min, x_max] straddle at least one point of the interleaved x range.
    {
        // Let F(x) := \int_{x_min}^{x_max} pdf(t, y) dt (y being the sample obtained above).
        // We create a temporary cdf values range, starting with 0 = F(x_min), ending with F(x_max),
        // and with all values for interleaved x in between.
        std::vector<Scalar> cdfx;
        cdfx.reserve(x_max.m_first - x_min.m_first + 2);
        cdfx.push_back(0.);

        Scalar const cdf0y = band_cdf(x_min, hit_y.m_first);
        for (auto x = x_min.m_first + 1; x <= x_max.m_first; ++x)
        {
            Interpol<std::vector<Scalar>> const xx(x, 0.);
            Scalar const cdfxy = band_cdf(xx, hit_y.m_first);
            cdfx.push_back(cdfxy - cdf0y);
        }
        Scalar const cdf1y = band_cdf(x_max, hit_y.m_first);
        cdfx.push_back(cdf1y - cdf0y);

        Interpol<std::vector<Scalar>> const xax = find_range(sx * cdfx.back(), cdfx);

        // Convert the interpol to one on the cdf positions vector.
        if(xax.m_first == cdfx.begin()) // Hit below first x-range position.
        {
            hit_x.m_first = x_min.m_first;
            hit_x.m_theta = x_min.m_theta + xax.m_theta * (1 - x_min.m_theta);
        }
        else if(xax.m_first == cdfx.end() - 2)
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

    result[0] = hit_x.value();

    pdf = get_pdf(hit_x, hit_y) / domain_proba;

    return result;
}

template<typename Scalar>
Scalar CDF2SamplerNonUniform<Scalar>::compute_pdf(vec2<Scalar> const& result,
        box2<Scalar> const& domain) const
{
    Interpol<std::vector<Scalar>> const x_min = find_range(domain[0][0], m_interleaved_x);
    Interpol<std::vector<Scalar>> const x_max = find_range(domain[1][0], m_interleaved_x);
    Interpol<std::vector<Scalar>> const y_min = find_range(domain[0][1], m_interleaved_y);
    Interpol<std::vector<Scalar>> const y_max = find_range(domain[1][1], m_interleaved_y);

    Scalar const cdf00 = get_cdf(x_min, y_min);
    Scalar const cdf10 = get_cdf(x_max, y_min);
    Scalar const cdf01 = get_cdf(x_min, y_max);
    Scalar const cdf11 = get_cdf(x_max, y_max);

    Scalar const domain_proba = (cdf11 - cdf10) + (cdf00 - cdf01);
    if(domain_proba < std::numeric_limits<Scalar>::epsilon())
    {
        return 0.;
    }

    Interpol<std::vector<Scalar>> const hit_x = find_range(result[0], m_interleaved_x);
    Interpol<std::vector<Scalar>> const hit_y = find_range(result[1], m_interleaved_y);

    return get_pdf(hit_x, hit_y) / domain_proba;
}

template class CDF2SamplerUniform<float> ;
template class CDF2SamplerUniform<double> ;
template class CDF2SamplerNonUniform<float> ;
template class CDF2SamplerNonUniform<double> ;

}
