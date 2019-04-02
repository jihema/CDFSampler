/*
 * \file cdf_tester.cpp
 */

#include "cdf_tester.h"
#include "cdf1_sampler.h"
#include "cdf2_sampler.h"
#include "random_stuff.h"
#include "scoped_timer.h"

#include <iostream>
#include <random>
#include <vector>

namespace cdf_sampler
{

namespace test
{

template<>
std::string CDFTester<CDF1SamplerUniform<float>>::get_sampler_name()
{
    return "CDF1SamplerUniform<float>";
}

template<>
std::string CDFTester<CDF1SamplerUniform<double>>::get_sampler_name()
{
    return "CDF1SamplerUniform<double>";
}

template<>
std::string CDFTester<CDF1SamplerNonUniform<float>>::get_sampler_name()
{
    return "CDF1SamplerNonUniform<float>";
}

template<>
std::string CDFTester<CDF1SamplerNonUniform<double>>::get_sampler_name()
{
    return "CDF1SamplerNonUniform<double>";
}

template<>
std::string CDFTester<CDF2SamplerUniform<float>>::get_sampler_name()
{
    return "CDF2SamplerUniform<float>";
}

template<>
std::string CDFTester<CDF2SamplerNonUniform<float>>::get_sampler_name()
{
    return "CDF2SamplerNonUniform<float>";
}

template<>
CDF1SamplerUniform<float> CDFTester<CDF1SamplerUniform<float>>::get_random_sampler()
{
    std::mt19937 random(m_seed);
    std::uniform_real_distribution<> uniform_01;

    long constexpr half_num_x = 2;
    std::vector<Scalar> f(2 * half_num_x + 1);

    for (Scalar& v : f)
    {
        v = 1. + 0.1 * uniform_01(random);
    }

    vec2<Scalar> domain { -M_PI_2, M_PI_2 };
    CDF1SamplerUniform<Scalar> sampler;
    sampler.init(domain, f);

    return sampler;
}

template<>
CDF1SamplerUniform<double> CDFTester<CDF1SamplerUniform<double>>::get_random_sampler()
{
    std::mt19937 random(m_seed);
    std::uniform_real_distribution<> uniform_01;

    long constexpr half_num_x = 2;
    std::vector<Scalar> f(2 * half_num_x + 1);

    for (Scalar& v : f)
    {
        v = 1. + 0.1 * uniform_01(random);
    }

    vec2<Scalar> domain { -M_PI_2, M_PI_2 };
    CDF1SamplerUniform<Scalar> sampler;
    sampler.init(domain, f);

    return sampler;
}

template<>
CDF1SamplerNonUniform<float> CDFTester<CDF1SamplerNonUniform<float>>::get_random_sampler()
{
    std::mt19937 random(m_seed);
    std::uniform_real_distribution<> uniform_01;

    long constexpr half_num_x = 2;
    std::vector<Scalar> f(2 * half_num_x + 1);

    for (Scalar& v : f)
    {
        v = 1. + 0.1 * uniform_01(random);
    }

    vec2<Scalar> domain { -M_PI_2, M_PI_2 };
    std::vector<Scalar> x = wiggle(domain[0], domain[1], 2 * half_num_x + 1, random);

    CDF1SamplerNonUniform<Scalar> sampler;
    sampler.init(x, f);

    return sampler;
}

template<>
CDF1SamplerNonUniform<double> CDFTester<CDF1SamplerNonUniform<double>>::get_random_sampler()
{
    std::mt19937 random(m_seed);
    std::uniform_real_distribution<> uniform_01;

    long constexpr half_num_x = 2;
    std::vector<Scalar> f(2 * half_num_x + 1);

    for (Scalar& v : f)
    {
        v = 1. + 0.1 * uniform_01(random);
    }

    vec2<Scalar> domain { -M_PI_2, M_PI_2 };
    std::vector<Scalar> x = wiggle(domain[0], domain[1], 2 * half_num_x + 1, random);

    CDF1SamplerNonUniform<Scalar> sampler;
    sampler.init(x, f);

    return sampler;
}

template<>
CDF2SamplerUniform<float> CDFTester<CDF2SamplerUniform<float>>::get_random_sampler()
{
    std::mt19937 random(m_seed);
    std::uniform_real_distribution<> uniform_01;

    size_t const x_size = 152;
    size_t const y_size = 341;
    std::vector<Scalar> f(x_size * y_size);
    for (Scalar& v : f)
    {
        v = 1. + 1. * uniform_01(random);
    }

    box2<Scalar> const domain { -2., -2., 2., 2. };
    CDF2SamplerUniform<Scalar> sampler;
    sampler.init(domain, x_size, y_size, f);

    return sampler;
}

template<>
CDF2SamplerNonUniform<float> CDFTester<CDF2SamplerNonUniform<float>>::get_random_sampler()
{
    std::mt19937 random(m_seed);
    std::uniform_real_distribution<> uniform_01;

    size_t const x_size = 152;
    size_t const y_size = 341;
    std::vector<Scalar> f(x_size * y_size);
    for (Scalar& v : f)
    {
        v = 1. + 1. * uniform_01(random);
    }

    box2<Scalar> const domain { -2., -2., 2., 2. };

    std::vector<Scalar> const x = wiggle(domain[0][0], domain[1][0], x_size, random);
    std::vector<Scalar> const y = wiggle(domain[0][1], domain[1][1], y_size, random);
    CDF2SamplerNonUniform<Scalar> sampler;
    sampler.init(x, y, f);

    return sampler;
}

template<typename Sampler>
bool CDFTester<Sampler>::test_random()
{
    std::cout << "Testing randomised " << get_sampler_name() << '\n';

    std::random_device seeder;  //Will be used to obtain a seed for the random number engine.
    std::random_device::result_type const seed = seeder();

    CDFTester tester(seed);
    Sampler const sampler = tester.get_random_sampler();
    typename Sampler::Domain const& domain = Sampler::max_domain;
    size_t const num_samples = 10000;

    tester.check_sum_one(sampler, domain, num_samples);

    std::cout << std::endl;

    return true;
}

template<typename Sampler>
void CDFTester<Sampler>::check_sum_one(Sampler const& sampler, typename Sampler::Domain const& domain,
        size_t const num_samples)
{
    ScopedTimer timer([](long ns)
    {   std::cout << "Time        = " << ns * 1.e-6 << " ms" << '\n';});

    using Domain = typename Sampler::Domain;
    using Scalar = typename Sampler::Scalar;

    std::mt19937 random(m_seed); // Standard mersenne_twister_engine seeded with rd().

    Scalar sum = 0.;
    for (size_t i = 0; i < num_samples; ++i)
    {
        auto const in_s = RandomVector<Scalar, Sampler::dimension>::get(random);
        Scalar pdf;
        sampler.sample(in_s, pdf, domain);
        sum += 1. / pdf;
    }

    Domain const& definition_domain = sampler.definition_domain();
    Scalar const domain_size = compute_area<Domain, Scalar>(intersect(domain, definition_domain));

    std::cout << "# samples   = " << num_samples << '\n';
    std::cout << "MC integral = " << sum / num_samples << '\n';
    std::cout << "Reference   = " << domain_size << '\n';
    std::cout << "Rel. error  = " << std::fabs(domain_size - sum / num_samples) / domain_size << '\n';
}

//template<typename Sampler>
//void CDFTester<Sampler>::check_sum_sin(Sampler const& sampler, typename Sampler::Domain const& domain,
//        size_t const num_samples)
//{
//    ScopedTimer timer([](long ns)
//    {   std::cout << "Time        = " << ns * 1.e-6 << " ms" << '\n';});
//
//    using Domain = typename Sampler::Domain;
//    using Scalar = typename Sampler::Scalar;
//
//    std::mt19937 random(m_seed); // Standard mersenne_twister_engine seeded with rd().
//
//    Scalar sum = 0.;
//    Scalar abs_sum = 0.;
//    for (size_t i = 0; i < num_samples; ++i)
//    {
//        auto const in_s = RandomVector<Scalar, 2>::get(random);
//        Scalar pdf;
//        vec2<Scalar> const xy = sampler.sample(in_s, pdf, domain);
//        sum += sin(xy[0]) * sin(xy[1]) / pdf;
//        abs_sum += fabs(sin(xy[0]) * sin(xy[1])) / pdf;
//    }
//
//    Scalar const reference = (cos(domain[0][0]) - cos(domain[1][0])) * (cos(domain[0][1]) - cos(domain[1][1]));
//
//    std::cout << "# samples   = " << num_samples << '\n';
//    std::cout << "MC integral = " << sum / num_samples << '\n';
//    std::cout << "Reference   = " << reference << '\n';
//    std::cout << "Rel. error  = " << std::fabs(reference - sum / num_samples) / (abs_sum / num_samples) << '\n';
//}

template class CDFTester<CDF1SamplerNonUniform<double>> ;
template class CDFTester<CDF1SamplerNonUniform<float>> ;
template class CDFTester<CDF1SamplerUniform<double>> ;
template class CDFTester<CDF1SamplerUniform<float>> ;
template class CDFTester<CDF2SamplerUniform<float>> ;
template class CDFTester<CDF2SamplerNonUniform<float>> ;

}

}
