/**
 * \file main.cpp
 */

#include "cdf_tester.h"
#include "cdf1_sampler.h"
#include "cdf2_sampler.h"

using namespace cdf_sampler;
using namespace cdf_sampler::test;

int main(int argc, char* argv[])
{
    CDFTester<CDF1SamplerUniform<float>>::test_random();
    CDFTester<CDF1SamplerNonUniform<float>>::test_random();
    CDFTester<CDF1SamplerUniform<double>>::test_random();
    CDFTester<CDF1SamplerNonUniform<double>>::test_random();
    CDFTester<CDF2SamplerUniform<float>>::test_random();
    CDFTester<CDF2SamplerNonUniform<float>>::test_random();
}
