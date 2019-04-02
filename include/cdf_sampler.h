/**
 * \file cdf_sampler.h
 */

#pragma once

#include <stddef.h>
#include <limits>
#include <vector>

namespace cdf_sampler
{

/**
 * \brief Base class for cdf samplers in any dimension.
 */
template<typename ScalarT>
class CDFSampler
{
public:

    using Scalar = ScalarT;

protected:

    CDFSampler() :
            m_sum(std::numeric_limits<Scalar>::signaling_NaN())
    {
    }

    virtual ~CDFSampler()
    {
    }

    /**
     * Computes an augmented vector of positions whose end points are the end points of x
     * and the (x.size() - 1) inner points are the middle of intervals in x.
     */
    static inline std::vector<Scalar> interleave(std::vector<Scalar> const& x)
    {
        std::vector<Scalar> interleaved(x.size() + 1);

        interleaved[0] = x[0];
        for (size_t ix = 1; ix < x.size(); ++ix)
        {
            interleaved[ix] = 0.5 * (x[ix - 1] + x[ix]);
        }
        interleaved.back() = x.back();

        return interleaved;
    }

    /**
     * After the cdf has been computed, this divides f and F by the integral, which is kept in m_sum.
     */
    inline void normalize()
    {
        m_sum = m_cdf.back();

        for (size_t i = 0; i < m_pdf.size(); ++i)
        {
            m_pdf[i] /= m_sum;
            m_cdf[i] /= m_sum;
        }
    }

    /**
     * Contains the pdf, integral normalised to 1.
     */
    std::vector<Scalar> m_pdf;

    /**
     * Contains the cdf, back() normalised to 1.
     */
    std::vector<Scalar> m_cdf;

    /**
     * Records the integral of the original input function (a.k.a. pdf before normalisation).
     */
    Scalar m_sum;
};

}
