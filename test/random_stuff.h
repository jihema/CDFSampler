/**
 * \file random_stuff.h
 */

#pragma once

#include "vec.h"

#include <random>

namespace cdf_sampler
{

namespace test
{

template<typename Scalar, int n>
class RandomVector;

template<typename Scalar>
class RandomVector<Scalar, 1>
{
public:
    static Scalar get(std::mt19937& random)
    {
        std::uniform_real_distribution<> uniform_01;
        return uniform_01(random);
    }
};

template<typename Scalar>
class RandomVector<Scalar, 2>
{
public:
    static vec2<Scalar> get(std::mt19937& random)
    {
        std::uniform_real_distribution<Scalar> uniform_01;
        return vec2 < Scalar > {uniform_01(random), uniform_01(random)};
    }
};

template<typename Scalar>
class RandomVector<Scalar, 3>
{
public:
    static vec3<Scalar> get(std::mt19937& random)
    {
        std::uniform_real_distribution<Scalar> uniform_01;
        return vec3<Scalar> { uniform_01(random), uniform_01(random),
                uniform_01(random) };
    }
};

}

}
