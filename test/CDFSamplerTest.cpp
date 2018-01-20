/**
 * \file RectifiedSamplingTest.cpp
 */

#include "CDFSamplerTest.h"

#include <cmath>
#include <iostream>

#include "../src/cdf1_sampler.h"
#include "../src/cdf2_sampler.h"
#include "../src/vec.h"
#include "random_stuff.h"
#include "scoped_timer.h"

namespace cdf_sampler
{

namespace test
{

template<typename Scalar>
void CDFSamplerTest<Scalar>::test_cdf_sampler_1D()
{
    std::random_device rd; //Will be used to obtain a seed for the random number engine
    std::mt19937 random(rd()); //Standard mersenne_twister_engine seeded with rd()
    std::uniform_real_distribution<> uniform_01;

    long constexpr half_num_x = 2;
    std::vector<Scalar> x(2 * half_num_x + 1);
    std::vector<Scalar> f(2 * half_num_x + 1);

    for (long i = -half_num_x; i <= half_num_x; ++i)
    {
        x[i + half_num_x] = i * M_PI_2 / half_num_x;
//        f[i + half_num_x] = 1 + std::cos(x[i + half_num_x]);
    }

    for (Scalar& v : f)
    {
        v = 1. + 0.1 * uniform_01(random);
    }

    vec2<Scalar> const& domain = CDF1Sampler<Scalar>::max_domain;
    size_t const num_samples = 100000;

    CDF1SamplerUniform<Scalar> cdf_sampler_uniform;
    cdf_sampler_uniform.init(vec2<Scalar> { -M_PI_2, M_PI_2 }, f);
    std::cout << "Uniform:\n";
    check_sum(cdf_sampler_uniform, domain, num_samples);

    CDF1SamplerNonUniform<Scalar> cdf_sampler_non_uniform;
    cdf_sampler_non_uniform.init(x, f);
    std::cout << "Non uniform:\n";
    check_sum(cdf_sampler_non_uniform, domain, num_samples);
}

template<typename Scalar>
void CDFSamplerTest<Scalar>::test_cdf_sampler_1D_again()
{
    std::random_device rd; //Will be used to obtain a seed for the random number engine
    std::mt19937 random(rd()); //Standard mersenne_twister_engine seeded with rd()
    std::uniform_real_distribution<> uniform_01;

    CDF1SamplerUniform<Scalar> cdf_sampler;
    long constexpr half_num_x = 4;
    std::vector<Scalar> x(2 * half_num_x + 1);
    std::vector<Scalar> f(2 * half_num_x + 1);
    for (long i = -half_num_x; i <= half_num_x; ++i)
    {
        x[i + half_num_x] = i * M_PI_2 / half_num_x;
        f[i + half_num_x] = 1;
    }
    cdf_sampler.init(vec2<Scalar> { -M_PI_2, M_PI_2 }, f);

    size_t const num_j = 100;
    Scalar const step = .2 / num_j;
    for (size_t j = 0; j <= num_j; ++j)
    {
        Scalar const x = .5 + j * step;
        Scalar const width = 0.05;
        vec2<Scalar> const domain { x, x + width };

        Scalar sum = 0.;
        size_t constexpr num_samples = 10000;
        for (size_t i = 0; i < num_samples; ++i)
        {
            Scalar const in_s = uniform_01(random);
            Scalar pdf;
            Scalar const out_s = cdf_sampler.sample(in_s, pdf, domain);
            sum += out_s / pdf;
        }

        std::cout << "{" << x << ", " << sum / num_samples << "}," << '\n';
    }
}

template<typename Scalar>
void CDFSamplerTest<Scalar>::test_cdf_sampler_2D()
{
    std::random_device rd; //Will be used to obtain a seed for the random number engine
    std::mt19937 random(rd()); //Standard mersenne_twister_engine seeded with rd()
    std::uniform_real_distribution<Scalar> uniform_01;

    size_t const x_size = 152;
    size_t const y_size = 341;
    std::vector<Scalar> f(x_size * y_size);
    for (Scalar& v : f)
    {
        v = 1. + 1. * uniform_01(random);
    }

    box2<Scalar> const definition_domain { -2., -2., 2., 2. };
    box2<Scalar> const sampling_domain { -uniform_01(random), -uniform_01(
            random), uniform_01(random), uniform_01(random) };
    size_t const num_samples = 1000000;
    unsigned int const seed = rd();

    {
        CDF2SamplerUniform<Scalar> cdf_sampler_uniform;
        cdf_sampler_uniform.init(definition_domain, x_size, y_size, f);
        std::cout << "Uniform sampler:\n";
        check_sum(cdf_sampler_uniform, sampling_domain, num_samples, &seed);
    }
    std::cout << '\n';
    {
        CDF2SamplerNonUniform<Scalar> cdf_sampler_non_uniform;
        cdf_sampler_non_uniform.init(definition_domain, x_size, y_size, f);
        std::cout << "Non-uniform sampler on uniform grid:\n";
        check_sum(cdf_sampler_non_uniform, sampling_domain, num_samples, &seed);
    }
    std::cout << '\n';
    {
        std::vector<Scalar> const x = wiggle(definition_domain[0][0],
                definition_domain[1][0], x_size, random);
        std::vector<Scalar> const y = wiggle(definition_domain[0][1],
                definition_domain[1][1], y_size, random);
        CDF2SamplerNonUniform<Scalar> cdf_sampler_non_uniform;
        cdf_sampler_non_uniform.init(x, y, f);
        std::cout << "Non-uniform sampler on non-uniform grid:\n";
        check_sum(cdf_sampler_non_uniform, sampling_domain, num_samples, &seed);
    }
    std::cout << '\n';
    {
        CDF2SamplerUniform<Scalar> uniform_sampler;
        uniform_sampler.init(definition_domain, x_size, y_size, f);
        std::cout << "Uniform sampler sin(x)*sin(y):\n";

        ScopedTimer timer([](long ns)
        {   std::cout << "Time        = " << ns * 1.e-6 << " ms" << '\n';});

        std::random_device rd; // Will be used to obtain a seed for the random number engine.
        std::mt19937 random(rd()); // Standard mersenne_twister_engine seeded with rd().

        Scalar sum = 0.;
        Scalar abs_sum = 0.;
        for (size_t i = 0; i < num_samples; ++i)
        {
            auto const in_s = RandomVector<Scalar, 2>::get(random);
            Scalar pdf;
            vec2<Scalar> const xy = uniform_sampler.sample(in_s, pdf,
                    sampling_domain);
            sum += sin(xy[0]) * sin(xy[1]) / pdf;
            abs_sum += fabs(sin(xy[0]) * sin(xy[1])) / pdf;
        }

        Scalar const reference = (cos(sampling_domain[0][0])
                - cos(sampling_domain[1][0]))
                * (cos(sampling_domain[0][1]) - cos(sampling_domain[1][1]));
        std::cout << "Reference   = " << reference << '\n';
        std::cout << "MC integral = " << sum / num_samples << '\n';
        std::cout << "Rel. error  = "
                << std::fabs(reference - sum / num_samples)
                        / (abs_sum / num_samples) << '\n';
        std::cout << "1/sqrt(n)   = " << 1. / std::sqrt(num_samples) << '\n';
    }
}

template<typename Scalar>
template<typename Sampler>
void CDFSamplerTest<Scalar>::check_sum(Sampler const& sampler,
        typename Sampler::Domain const& domain, size_t const num_samples,
        unsigned int const* const seed)
{
    ScopedTimer timer([](long ns)
    {   std::cout << "Time        = " << ns * 1.e-6 << " ms" << '\n';});

    using Domain = typename Sampler::Domain;

    std::random_device rd; // Will be used to obtain a seed for the random number engine.
    std::mt19937 random(seed ? *seed : rd()); // Standard mersenne_twister_engine seeded with rd().

    Scalar sum = 0.;
    for (size_t i = 0; i < num_samples; ++i)
    {
        auto const in_s = RandomVector<Scalar, Sampler::dimension>::get(random);
        Scalar pdf;
        sampler.sample(in_s, pdf, domain);
        sum += 1. / pdf;
    }

    Domain const& definition_domain = sampler.definition_domain();
    Scalar const domain_size = compute_area<Domain, Scalar>(
            intersect(domain, definition_domain));

    std::cout << "Reference   = " << domain_size << '\n';
    std::cout << "MC integral = " << sum / num_samples << '\n';
    std::cout << "Rel. error  = "
            << std::fabs(domain_size - sum / num_samples) / domain_size << '\n';
    std::cout << "1/sqrt(n)   = " << 1. / std::sqrt(num_samples) << '\n';
}

template class CDFSamplerTest<float> ;
template class CDFSamplerTest<double> ;

}

}
