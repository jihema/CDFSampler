/**
 * \file cdf2_sampler.h
 */

#pragma once

#include <vector>
#include <gmath_bbox2.h>
#include <gmath_vec2.h>
#include <array>
#include <cmath>

#include "interpol.h"
#include "cdf_sampler.h"

namespace dneg
{

template<typename Scalar>
class CDF2Sampler: public CDFSampler<Scalar>
{
public:
    using Domain = GMathBbox2<Scalar>;
    static constexpr int dimension = 2;

    virtual void init(Domain const& xy_range, size_t const x_size, size_t const y_size,
            std::vector<Scalar> const& f) = 0;

    virtual GMathVec2<Scalar> sample(GMathVec2<Scalar> const& sxy, Scalar& pdf,
            GMathBbox2<Scalar> const& domain) const = 0;

// FIXME this should work with protected:
    inline Scalar& get_cdf(size_t xi, size_t yi)
    {
        return this->m_cdf[yi * m_x_size + xi];
    }

    inline Scalar get_cdf(size_t xi, size_t yi) const
    {
        return this->m_cdf[yi * m_x_size + xi];
    }

    inline Scalar get_pdf(size_t xi, size_t yi) const
    {
        return this->m_pdf[yi * m_x_size + xi];
    }

    virtual GMathVec2<Scalar> base_xy(size_t xi, size_t yi) const = 0;
    virtual GMathVec2<Scalar> interleaved_xy(size_t xi, size_t yi) const = 0;

    void print_cdf(std::ostream& os) const
    {
        {
            GMathVec2<Scalar> const xy = interleaved_xy(0, 0);
            os << "{" << xy[0] << ", " << xy[1] << ", " << 0. << "},\n";
        }

        for (size_t xi = 0; xi < m_x_size; ++xi)
        {
            GMathVec2<Scalar> const xy = interleaved_xy(xi + 1, 0);
            os << "{" << xy[0] << ", " << xy[1] << ", " << 0. << "},\n";
        }

        for (size_t yi = 0; yi < m_y_size; ++yi)
        {
            GMathVec2<Scalar> const xy = interleaved_xy(0, yi + 1);
            os << "{" << xy[0] << ", " << xy[1] << ", " << 0. << "},\n";
        }

        for (size_t yi = 0; yi < m_y_size; ++yi)
        {
            for (size_t xi = 0; xi < m_x_size; ++xi)
            {
                GMathVec2<Scalar> const xy = interleaved_xy(xi + 1, yi + 1);
                os << "{" << xy[0] << ", " << xy[1] << ", " << get_cdf(xi, yi) << "},\n";
            }
        }
    }

    void print_pdf(std::ostream& os) const
    {
        for (size_t yi = 0; yi < m_y_size; ++yi)
        {
            for (size_t xi = 0; xi < m_x_size; ++xi)
            {
                GMathVec2<Scalar> const xy = base_xy(xi, yi);
                os << "{" << xy[0] << ", " << xy[1] << ", " << get_pdf(xi, yi) << ",\n";
            }
        }
    }

protected:

    inline std::array<Scalar, 4> corner_cdfs(size_t const xi, size_t const yi) const
    {
        std::array < Scalar, 4 > cdfs;

        // In xy order: 00, 10, 01, 11.
        cdfs[0] = (xi == 0 || yi == 0) ? 0. : get_cdf(xi - 1, yi - 1);
        cdfs[1] = yi == 0 ? 0. : get_cdf(xi, yi - 1);
        cdfs[2] = xi == 0 ? 0. : get_cdf(xi - 1, yi);
        cdfs[3] = get_cdf(xi, yi);

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

    GMathVec2<Scalar> sample(GMathVec2<Scalar> const& sxy, Scalar& pdf,
            GMathBbox2<Scalar> const& domain) const override;

    Scalar interleaved_x(size_t xi) const
    {
        Scalar x;
        if(xi == 0)
        {
            x = m_x_min;
        }
        else if(xi == this->m_x_size)
        {
            x = m_x_max;
        }
        else
        {
            x = m_x_min + (xi - 0.5) * (m_x_max - m_x_min) / (this->m_x_size - 1);
        }

        return x;
    }

    Scalar interleaved_y(size_t yi) const
    {
        Scalar y;
        if(yi == 0)
        {
            y = m_y_min;
        }
        else if(yi == this->m_y_size)
        {
            y = m_y_max;
        }
        else
        {
            y = m_y_min + (yi - 0.5) * (m_y_max - m_y_min) / (this->m_y_size - 1);
        }

        return y;
    }

    inline GMathVec2<Scalar> interleaved_xy(size_t xi, size_t yi) const override
    {
        Scalar x;
        if(xi == 0)
        {
            x = m_x_min;
        }
        else if(xi == this->m_x_size)
        {
            x = m_x_max;
        }
        else
        {
            x = m_x_min + (xi - 0.5) * (m_x_max - m_x_min) / (this->m_x_size - 1);
        }

        Scalar y;
        if(yi == 0)
        {
            y = m_y_min;
        }
        else if(yi == this->m_y_size)
        {
            y = m_y_max;
        }
        else
        {
            y = m_y_min + (yi - 0.5) * (m_y_max - m_y_min) / (this->m_y_size - 1);
        }

        return GMathVec2<Scalar>(x, y);
    }

    GMathVec2<Scalar> base_xy(size_t xi, size_t yi) const
    {
        Scalar const x = m_x_min + xi * (m_x_max - m_x_min) / (this->m_x_size - 1);
        Scalar const y = m_y_min + yi * (m_y_max - m_y_min) / (this->m_y_size - 1);
        return GMathVec2<Scalar>(x, y);
    }

    Domain definition_domain() const
    {
        return Domain(m_x_min, m_y_min, m_x_max, m_y_max);
    }

private:

    inline Scalar get_cdf(size_t xi, size_t yi) const
    {
        return this->m_cdf[yi * this->m_x_size + xi];
    }

    inline Scalar get_pdf(size_t xi, size_t yi) const
    {
        return this->m_pdf[yi * this->m_x_size + xi];
    }

    inline Scalar band_cdf(size_t const ix, Scalar const x_theta, Scalar const iy) const
    {
        std::array<Scalar, 4> const cdfs = this->corner_cdfs(ix, iy); // In xy order: 00, 10, 01, 11.
        return linear_interpolation(cdfs[2] - cdfs[0], cdfs[3] - cdfs[1], x_theta);
    }

    inline Scalar get_cdf(Scalar const x, Scalar const y, size_t& ix, Scalar& theta_x, size_t& iy,
            Scalar& theta_y) const
    {
        assert(x >= m_x_min && x <= m_x_max);
        size_t const x_size = this->m_x_size;

        if(x == m_x_max)
        {
            ix = x_size - 1;
            theta_x = 1.;
        }

        Scalar const scaled_x = (x_size - 1) * (x - m_x_min) / (m_x_max - m_x_min); // in [0, m_pdf.size()-1).

        ix = std::floor(scaled_x + 0.5); // in {0, ..., m_cdf.size()-1}

        if(ix == 0)
        {
            theta_x = 2. * scaled_x;
        }
        else if(ix == x_size - 1)
        {
            theta_x = 2 * (scaled_x - ix) + 1.;
        }
        else
        {
            theta_x = scaled_x - ix + 0.5;
        }

        assert(y >= m_y_min && y <= m_y_max);
        size_t const y_size = this->m_y_size;

        if(y == m_y_max)
        {
            iy = y_size - 1;
            theta_y = 1.;
        }

        Scalar const scaled_y = (y_size - 1) * (y - m_y_min) / (m_y_max - m_y_min); // in [0, m_pdf.size()-1).

        iy = std::floor(scaled_y + 0.5); // in {0, ..., m_cdf.size()-1}

        if(iy == 0)
        {
            theta_y = 2. * scaled_y;
        }
        else if(iy == y_size - 1)
        {
            theta_y = 2 * (scaled_y - iy) + 1.;
        }
        else
        {
            theta_y = scaled_y - iy + 0.5;
        }

        std::array < Scalar, 4 > cdfs = this->corner_cdfs(ix, iy);

        return bilinear_interpolation(cdfs[0], cdfs[1], cdfs[2], cdfs[3], theta_x, theta_y);
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
     * @param x Grid projections on the x-axis.
     * @param y Grid projections on the y-axis.
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
            Scalar const step = (xy_range[1][0] - xy_range[0][0]) / (x_size - 1);
            for (size_t i = 0; i < x_size; ++i)
            {
                x[i] = xy_range[0][0] + i * step;
            }
        }

        std::vector<Scalar> y(y_size);
        {
            Scalar const step = (xy_range[1][1] - xy_range[0][1]) / (y_size - 1);
            for (size_t i = 0; i < y_size; ++i)
            {
                y[i] = xy_range[0][1] + i * step;
            }
        }

        init(x, y, f);
    }

    GMathVec2<Scalar> sample(GMathVec2<Scalar> const& sxy, Scalar& pdf,
            GMathBbox2<Scalar> const& domain) const override;

    Scalar compute_pdf(GMathVec2<Scalar> const& result,
            GMathBbox2<Scalar> const& rectified_domain) const;

    GMathVec2<Scalar> base_xy(size_t xi, size_t yi) const override
    {
        return GMathVec2<Scalar>(m_x[xi], m_y[yi]);
    }

    GMathVec2<Scalar> interleaved_xy(size_t xi, size_t yi) const override
    {
        return GMathVec2<Scalar>(m_interleaved_x[xi], m_interleaved_y[yi]);
    }

    Domain definition_domain() const
    {
        return Domain(m_x.front(), m_y.front(), m_x.back(), m_y.back());
    }

private:

    inline Scalar band_cdf(Interpol<std::vector<Scalar>> const& x,
            typename std::vector<Scalar>::const_iterator const& y) const
    {
        assert(x.m_first >= m_interleaved_x.begin() && x.m_first < m_interleaved_x.end() - 1);
        assert(y >= m_interleaved_y.begin() && y < m_interleaved_y.end() - 1);

        size_t const xi = x.m_first - m_interleaved_x.begin();
        size_t const yi = y - m_interleaved_y.begin();

        std::array<Scalar, 4> const cdfs = this->corner_cdfs(xi, yi); // In xy order: 00, 10, 01, 11.
        return linear_interpolation(cdfs[2] - cdfs[0], cdfs[3] - cdfs[1], x.m_theta);
    }

    inline Scalar get_cdf(Interpol<std::vector<Scalar>> const& x,
            Interpol<std::vector<Scalar>> const& y) const
    {
        assert(x.m_first >= m_interleaved_x.begin() && x.m_first < m_interleaved_x.end() - 1);
        assert(y.m_first >= m_interleaved_y.begin() && y.m_first < m_interleaved_y.end() - 1);

        size_t const xi = x.m_first - m_interleaved_x.begin();
        size_t const yi = y.m_first - m_interleaved_y.begin();

        std::array<Scalar, 4> const cdfs = this->corner_cdfs(xi, yi); // In xy order: 00, 10, 01, 11.
        return bilinear_interpolation(cdfs[0], cdfs[1], cdfs[2], cdfs[3], x.m_theta, y.m_theta);
    }

    inline Scalar get_pdf(Interpol<std::vector<Scalar>> const& hit_x,
            Interpol<std::vector<Scalar>> const& hit_y) const
    {
        assert(m_interleaved_x.begin() <= hit_x.m_first && hit_x.m_first < m_interleaved_x.end());
        assert(m_interleaved_y.begin() <= hit_y.m_first && hit_y.m_first < m_interleaved_y.end());

        size_t const idx_x = hit_x.m_first - m_interleaved_x.begin();
        size_t const idx_y = hit_y.m_first - m_interleaved_y.begin();

        return static_cast<CDF2Sampler<Scalar> const* const >(this)->get_pdf(idx_x, idx_y);
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
