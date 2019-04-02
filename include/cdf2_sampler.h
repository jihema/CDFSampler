/**
 * \file cdf2_sampler.h
 */

#pragma once

#include <stddef.h>
#include <array>
#include <cassert>
#include <iostream>
#include <vector>
#include <cmath>

#include "cdf_sampler.h"
#include "interpol.h"
#include "vec.h"

namespace cdf_sampler
{

template<typename Scalar>
class CDF2Sampler: public CDFSampler<Scalar>
{
public:
    virtual ~CDF2Sampler()
    {
    }

    using Domain = box2<Scalar>;
    static constexpr int dimension = 2;
    static Domain const max_domain;

    virtual void init(Domain const& xy_range, size_t const x_size,
            size_t const y_size, std::vector<Scalar> const& f) = 0;

    virtual vec2<Scalar> sample(vec2<Scalar> const& in_s, Scalar& pdf,
            box2<Scalar> const& domain) const = 0;

// FIXME this should work with protected:
    inline Scalar& get_cdf(size_t ix, size_t iy)
    {
        return this->m_cdf[iy * m_x_size + ix];
    }

    inline Scalar get_cdf(size_t ix, size_t iy) const
    {
        return this->m_cdf[iy * m_x_size + ix];
    }

    inline Scalar get_pdf(size_t ix, size_t iy) const
    {
        return this->m_pdf[iy * m_x_size + ix];
    }

    virtual vec2<Scalar> base_xy(size_t ix, size_t iy) const = 0;
    virtual vec2<Scalar> interleaved_xy(size_t ix, size_t iy) const = 0;

    void print_cdf(std::ostream& os) const
    {
        {
            vec2<Scalar> const xy = interleaved_xy(0, 0);
            os << "{" << xy[0] << ", " << xy[1] << ", " << 0. << "},\n";
        }

        for (size_t ix = 0; ix < m_x_size; ++ix)
        {
            vec2<Scalar> const xy = interleaved_xy(ix + 1, 0);
            os << "{" << xy[0] << ", " << xy[1] << ", " << 0. << "},\n";
        }

        for (size_t iy = 0; iy < m_y_size; ++iy)
        {
            vec2<Scalar> const xy = interleaved_xy(0, iy + 1);
            os << "{" << xy[0] << ", " << xy[1] << ", " << 0. << "},\n";
        }

        for (size_t iy = 0; iy < m_y_size; ++iy)
        {
            for (size_t ix = 0; ix < m_x_size; ++ix)
            {
                vec2<Scalar> const xy = interleaved_xy(ix + 1, iy + 1);
                os << "{" << xy[0] << ", " << xy[1] << ", " << get_cdf(ix, iy)
                        << "},\n";
            }
        }
    }

    void print_pdf(std::ostream& os) const
    {
        for (size_t iy = 0; iy < m_y_size; ++iy)
        {
            for (size_t ix = 0; ix < m_x_size; ++ix)
            {
                vec2<Scalar> const xy = base_xy(ix, iy);
                os << "{" << xy[0] << ", " << xy[1] << ", " << get_pdf(ix, iy)
                        << ",\n";
            }
        }
    }

protected:

    inline std::array<Scalar, 4> corner_cdfs(size_t const ix,
            size_t const iy) const
    {
        std::array<Scalar, 4> cdfs;

        // In xy order: 00, 10, 01, 11.
        cdfs[0] = (ix == 0 || iy == 0) ? 0. : get_cdf(ix - 1, iy - 1);
        cdfs[1] = iy == 0 ? 0. : get_cdf(ix, iy - 1);
        cdfs[2] = ix == 0 ? 0. : get_cdf(ix - 1, iy);
        cdfs[3] = get_cdf(ix, iy);

        return cdfs;
    }

    size_t m_x_size;
    size_t m_y_size;
};

template<typename Scalar>
class CDF2SamplerUniform: public CDF2Sampler<Scalar>
{
public:
    using Domain = typename CDF2Sampler<Scalar>::Domain;

    virtual ~CDF2SamplerUniform()
    {
    }

    void init(Domain const& xy_range, size_t const x_size, size_t const y_size,
            std::vector<Scalar> const& f) override;

    vec2<Scalar> sample(vec2<Scalar> const& sxy, Scalar& pdf,
            box2<Scalar> const& domain) const override;

    Scalar interleaved_x(size_t ix) const
    {
        Scalar x;
        if (ix == 0)
        {
            x = m_x_min;
        }
        else if (ix == this->m_x_size)
        {
            x = m_x_max;
        }
        else
        {
            x = m_x_min
                    + (ix - 0.5) * (m_x_max - m_x_min) / (this->m_x_size - 1);
        }

        return x;
    }

    Scalar interleaved_y(size_t iy) const
    {
        Scalar y;
        if (iy == 0)
        {
            y = m_y_min;
        }
        else if (iy == this->m_y_size)
        {
            y = m_y_max;
        }
        else
        {
            y = m_y_min
                    + (iy - 0.5) * (m_y_max - m_y_min) / (this->m_y_size - 1);
        }

        return y;
    }

    inline vec2<Scalar> interleaved_xy(size_t ix, size_t iy) const override
    {
        Scalar x;
        if (ix == 0)
        {
            x = m_x_min;
        }
        else if (ix == this->m_x_size)
        {
            x = m_x_max;
        }
        else
        {
            x = m_x_min
                    + (ix - 0.5) * (m_x_max - m_x_min) / (this->m_x_size - 1);
        }

        Scalar y;
        if (iy == 0)
        {
            y = m_y_min;
        }
        else if (iy == this->m_y_size)
        {
            y = m_y_max;
        }
        else
        {
            y = m_y_min
                    + (iy - 0.5) * (m_y_max - m_y_min) / (this->m_y_size - 1);
        }

        return vec2<Scalar> { x, y };
    }

    vec2<Scalar> base_xy(size_t ix, size_t iy) const override
    {
        Scalar const x = m_x_min
                + ix * (m_x_max - m_x_min) / (this->m_x_size - 1);
        Scalar const y = m_y_min
                + iy * (m_y_max - m_y_min) / (this->m_y_size - 1);
        return vec2<Scalar> { x, y };
    }

    Domain definition_domain() const
    {
        return Domain { vec2<Scalar> { m_x_min, m_y_min }, vec2<Scalar> {
                m_x_max, m_y_max } };
    }

private:

    inline Scalar get_cdf(size_t ix, size_t iy) const
    {
        return this->m_cdf[iy * this->m_x_size + ix];
    }

    inline Scalar get_cdf(size_t const ix, Scalar const theta_x,
            size_t const iy, Scalar const theta_y) const
    {
        std::array<Scalar, 4> const cdfs = this->corner_cdfs(ix, iy); // In xy order: 00, 10, 01, 11.
        return bilinear_interpolation(cdfs[0], cdfs[1], cdfs[2], cdfs[3],
                theta_x, theta_y);
    }

    inline Scalar get_pdf(size_t ix, size_t iy) const
    {
        return this->m_pdf[iy * this->m_x_size + ix];
    }

    inline Scalar band_cdf(size_t const ix, Scalar const x_theta,
            Scalar const iy) const
    {
        std::array<Scalar, 4> const cdfs = this->corner_cdfs(ix, iy); // In xy order: 00, 10, 01, 11.
        return linear_interpolation(cdfs[2] - cdfs[0], cdfs[3] - cdfs[1],
                x_theta);
    }

    inline Scalar get_cdf(Scalar const x, Scalar const y, size_t& ix,
            Scalar& theta_x, size_t& iy, Scalar& theta_y) const
    {
        assert(x >= m_x_min && x <= m_x_max);
        size_t const x_size = this->m_x_size;

        if (x == m_x_max)
        {
            ix = x_size - 1;
            theta_x = 1.;
        }

        Scalar const scaled_x = (x_size - 1) * (x - m_x_min)
                / (m_x_max - m_x_min); // in [0, m_pdf.size()-1).

        ix = std::floor(scaled_x + 0.5); // in {0, ..., m_cdf.size()-1}

        if (ix == 0)
        {
            theta_x = 2. * scaled_x;
        }
        else if (ix == x_size - 1)
        {
            theta_x = 2 * (scaled_x - ix) + 1.;
        }
        else
        {
            theta_x = scaled_x - ix + 0.5;
        }

        assert(y >= m_y_min && y <= m_y_max);
        size_t const y_size = this->m_y_size;

        if (y == m_y_max)
        {
            iy = y_size - 1;
            theta_y = 1.;
        }

        Scalar const scaled_y = (y_size - 1) * (y - m_y_min)
                / (m_y_max - m_y_min); // in [0, m_pdf.size()-1).

        iy = std::floor(scaled_y + 0.5); // in {0, ..., m_cdf.size()-1}

        if (iy == 0)
        {
            theta_y = 2. * scaled_y;
        }
        else if (iy == y_size - 1)
        {
            theta_y = 2 * (scaled_y - iy) + 1.;
        }
        else
        {
            theta_y = scaled_y - iy + 0.5;
        }

        std::array<Scalar, 4> cdfs = this->corner_cdfs(ix, iy);

        return bilinear_interpolation(cdfs[0], cdfs[1], cdfs[2], cdfs[3],
                theta_x, theta_y);
    }

    Scalar m_x_min;
    Scalar m_y_min;
    Scalar m_x_max;
    Scalar m_y_max;
};

template<typename Scalar>
class CDF2SamplerNonUniform: public CDF2Sampler<Scalar>
{
public:
    using Domain = typename CDF2Sampler<Scalar>::Domain;

    virtual ~CDF2SamplerNonUniform()
    {
    }

    /**
     * \brief Builds a 2D cdf over a (non-uniform) Cartesian grid.
     *
     * @param x Grid projections on the x-aixs.
     * @param y Grid projections on the y-aixs.
     * @param f A non-uniform function sampled over the (x,y) grid.
     */
    void init(std::vector<Scalar> const& x, std::vector<Scalar> const& y,
            std::vector<Scalar> const& f);

    void init(Domain const& xy_range, size_t const x_size, size_t const y_size,
            std::vector<Scalar> const& f) override
    {
        assert(f.size() == x_size * y_size);

        std::vector<Scalar> x(x_size);
        {
            Scalar const step = (xy_range[1][0] - xy_range[0][0])
                    / (x_size - 1);
            for (size_t i = 0; i < x_size; ++i)
            {
                x[i] = xy_range[0][0] + i * step;
            }
        }

        std::vector<Scalar> y(y_size);
        {
            Scalar const step = (xy_range[1][1] - xy_range[0][1])
                    / (y_size - 1);
            for (size_t i = 0; i < y_size; ++i)
            {
                y[i] = xy_range[0][1] + i * step;
            }
        }

        init(x, y, f);
    }

    vec2<Scalar> sample(vec2<Scalar> const& sxy, Scalar& pdf,
            box2<Scalar> const& domain) const override;

    Scalar compute_pdf(vec2<Scalar> const& result,
            box2<Scalar> const& rectified_domain) const;

    vec2<Scalar> base_xy(size_t ix, size_t iy) const override
    {
        return vec2<Scalar> { m_x[ix], m_y[iy] };
    }

    vec2<Scalar> interleaved_xy(size_t ix, size_t iy) const override
    {
        return vec2<Scalar> { m_interleaved_x[ix], m_interleaved_y[iy] };
    }

    Domain definition_domain() const
    {
        return Domain { vec2<Scalar> { m_x.front(), m_y.front() },
                vec2<Scalar> { m_x.back(), m_y.back() } };
    }

private:

    inline Scalar band_cdf(Interpol<std::vector<Scalar>> const& x,
            typename std::vector<Scalar>::const_iterator const& y) const
    {
        assert(
                x.m_first >= m_interleaved_x.begin()
                        && x.m_first < m_interleaved_x.end() - 1);
        assert(y >= m_interleaved_y.begin() && y < m_interleaved_y.end() - 1);

        size_t const ix = x.m_first - m_interleaved_x.begin();
        size_t const iy = y - m_interleaved_y.begin();

        std::array<Scalar, 4> const cdfs = this->corner_cdfs(ix, iy); // In xy order: 00, 10, 01, 11.
        return linear_interpolation(cdfs[2] - cdfs[0], cdfs[3] - cdfs[1],
                x.m_theta);
    }

    inline Scalar get_cdf(Interpol<std::vector<Scalar>> const& x,
            Interpol<std::vector<Scalar>> const& y) const
    {
        assert(
                x.m_first >= m_interleaved_x.begin()
                        && x.m_first < m_interleaved_x.end() - 1);
        assert(
                y.m_first >= m_interleaved_y.begin()
                        && y.m_first < m_interleaved_y.end() - 1);

        size_t const ix = x.m_first - m_interleaved_x.begin();
        size_t const iy = y.m_first - m_interleaved_y.begin();

        std::array<Scalar, 4> const cdfs = this->corner_cdfs(ix, iy); // In xy order: 00, 10, 01, 11.
        return bilinear_interpolation(cdfs[0], cdfs[1], cdfs[2], cdfs[3],
                x.m_theta, y.m_theta);
    }

    inline Scalar get_pdf(Interpol<std::vector<Scalar>> const& hit_x,
            Interpol<std::vector<Scalar>> const& hit_y) const
    {
        assert(
                m_interleaved_x.begin() <= hit_x.m_first
                        && hit_x.m_first < m_interleaved_x.end());
        assert(
                m_interleaved_y.begin() <= hit_y.m_first
                        && hit_y.m_first < m_interleaved_y.end());

        size_t const idx_x = hit_x.m_first - m_interleaved_x.begin();
        size_t const idx_y = hit_y.m_first - m_interleaved_y.begin();

        return static_cast<CDF2Sampler<Scalar> const* const >(this)->get_pdf(
                idx_x, idx_y);
    }

    std::vector<Scalar> m_x; ///< Positions on the x axis.
    std::vector<Scalar> m_y; ///< Positions on the y axis.

    /**
     *  Interleaved x values: {x0, (x0 + x1)/2, ..., (x_{n-2} + x_{n-1})/2, x_{n-1}}.
     *  Note that m_interleaved_x.size() = m_x.size()+1.
     */
    std::vector<Scalar> m_interleaved_x;

    /**
     *  Interleaved y values: {y0, (y0 + y1)/2, ..., (y_{n-2} + y_{n-1})/2, y_{n-1}}.
     *  Note that m_interleaved_y.size() = m_y.size()+1.
     */
    std::vector<Scalar> m_interleaved_y;

};

}
