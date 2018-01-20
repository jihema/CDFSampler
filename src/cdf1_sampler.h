/**
 * \file cdf1_sampler.h
 */

#pragma once

#include <stddef.h>
#include <cassert>
#include <vector>
#include <cmath>

#include "cdf_sampler.h"
#include "interpol.h"
#include "vec.h"

namespace cdf_sampler
{

template<typename Scalar>
class CDF1Sampler: public CDFSampler<Scalar>
{
public:

    using Domain = vec2<Scalar>;
    static Domain const max_domain;

    static constexpr int dimension = 1;

    virtual ~CDF1Sampler()
    {
    }

    virtual void init(Domain const& x_range, std::vector<Scalar> const& f) = 0;

    virtual Scalar sample(Scalar const in_s, Scalar& pdf,
            Domain const& domain) const = 0;

    virtual Scalar interleaved_x(size_t ix) const = 0;

    virtual Domain definition_domain() const = 0;

protected:

    inline Scalar& get_cdf(size_t ix)
    {
        return this->m_cdf[ix];
    }

    inline Scalar get_cdf(size_t ix) const
    {
        return this->m_cdf[ix];
    }

    inline Scalar get_pdf(size_t ix) const
    {
        return this->m_pdf[ix];
    }

};

template<typename Scalar>
class CDF1SamplerUniform: public CDF1Sampler<Scalar>
{
public:
    using Domain = typename CDF1Sampler<Scalar>::Domain;

    virtual ~CDF1SamplerUniform()
    {
    }

    void init(Domain const& x_range, std::vector<Scalar> const& f) override;

    Scalar sample(Scalar const s0, Scalar& pdf, Domain const& domain =
            CDF1Sampler<Scalar>::max_domain) const override;

    inline Scalar interleaved_x(size_t ix) const override
    {
        if (ix == 0)
        {
            return m_x_min;
        }
        else if (ix == this->m_cdf.size())
        {
            return m_x_max;
        }
        else
        {
            return m_x_min
                    + (ix - 0.5) * (m_x_max - m_x_min)
                            / (this->m_cdf.size() - 1);
        }
    }

    Domain definition_domain() const override
    {
        return Domain { m_x_min, m_x_max };
    }

private:

    /**
     * \brief Computes linearly interpolated cdf, sets interpolation parameters.
     *
     * @param x in [x_min, x_max].
     *
     * @param ix gets the index (0 <= ix < cdf.size()) in the interleaved range,
     * such that x \in [interleaved_x(xi), interleaved_x(xi + 1)) (except when x == x_max,
     * then x = interleaved_x(xi + 1)).
     *
     * @param theta gets the interpolation parameter for x in the above interval, in [0,1)
     * (except when x == x_max, then theta = 1).
     *
     * @return cdf value at x.
     */
    inline Scalar get_cdf(Scalar const x, size_t& ix, Scalar& theta) const
    {
        assert(x >= m_x_min && x <= m_x_max);

        if (x == m_x_max)
        {
            ix = this->m_cdf.size() - 1;
            theta = 1.;
            return this->m_cdf.back();
        }

        Scalar const scaled_x = (this->m_cdf.size() - 1) * (x - m_x_min)
                / (m_x_max - m_x_min); // in [0, m_pdf.size()-1).

        ix = std::floor(scaled_x + 0.5); // in {0, ..., m_cdf.size()-1}

        if (ix == 0)
        {
            theta = 2. * scaled_x;
        }
        else if (ix == this->m_cdf.size() - 1)
        {
            theta = 2 * (scaled_x - ix) + 1.;
        }
        else
        {
            theta = scaled_x - ix + 0.5;
        }

        return linear_interpolation(ix == 0 ? 0 : this->m_cdf[ix - 1],
                this->m_cdf[ix], theta);
    }

    Scalar m_x_min;
    Scalar m_x_max;
};

template<typename Scalar>
class CDF1SamplerNonUniform: public CDF1Sampler<Scalar>
{
public:
    using Domain = typename CDF1Sampler<Scalar>::Domain;

    virtual ~CDF1SamplerNonUniform()
    {
    }

    void init(vec2<Scalar> const& x_range, std::vector<Scalar> const& f)
            override
    {
        std::vector<Scalar> x(f.size());
        Scalar const step = (x_range[1] - x_range[0]) / (f.size() - 1);
        for (size_t i = 0; i < f.size(); ++i)
        {
            x[i] = x_range[0] + i * step;
        }

        init(x, f);
    }

    void init(std::vector<Scalar> const& x, std::vector<Scalar> const& f);

    Scalar sample(Scalar const in_s, Scalar& pdf, vec2<Scalar> const& domain =
            CDF1Sampler<Scalar>::max_domain) const override;

    inline Scalar interleaved_x(size_t ix) const override
    {
        return m_interleaved_x[ix];
    }

    inline Domain definition_domain() const override
    {
        return vec2<Scalar> { m_x.front(), m_x.back() };
    }

private:

    inline Scalar get_cdf(Interpol<std::vector<Scalar>> const& x) const
    {
        assert(
                x.m_first >= m_interleaved_x.begin()
                        && x.m_first < m_interleaved_x.end() - 1);

        size_t const xi = x.m_first - m_interleaved_x.begin();

        Scalar const cdf0 =
                (xi == 0) ? 0. : CDF1Sampler<Scalar>::get_cdf(xi - 1);
        Scalar const cdf1 = CDF1Sampler<Scalar>::get_cdf(xi);

        return linear_interpolation(cdf0, cdf1, x.m_theta);
    }

    inline Scalar get_pdf(Interpol<std::vector<Scalar>> const& x) const
    {
        assert(
                m_interleaved_x.begin() <= x.m_first
                        && x.m_first < m_interleaved_x.end());

        size_t const idx_x = x.m_first - m_interleaved_x.begin();

        return this->m_pdf[idx_x];
    }

    std::vector<Scalar> m_x; ///< Positions on the x axis.

    /**
     * \brief Interleaved x values: {x0, (x0 + x1)/2, ..., (x_{n-2} + x_{n-1})/2, x_{n-1}}.
     *
     * Note that m_interleaved_x.size() = m_x.size() + 1.
     *
     * In our model, the pdf is constant on any interval between two consecutive interleaved x.
     */
    std::vector<Scalar> m_interleaved_x;

};

}
