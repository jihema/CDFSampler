/**
 * \file random_stuff.h
 */

#pragma once

#include <random>
#include <gmath_vec2.h>
#include <gmath_vec3.h>

namespace dneg
{

using namespace light;

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
    static GMathVec2<Scalar> get(std::mt19937& random)
    {
        std::uniform_real_distribution<> uniform_01;
        return GMathVec2<Scalar>(uniform_01(random), uniform_01(random));
    }
};

template<typename Scalar>
class RandomVector<Scalar, 3>
{
public:
    static GMathVec3<Scalar> get(std::mt19937& random)
    {
        std::uniform_real_distribution<> uniform_01;
        return GMathVec3<Scalar>(uniform_01(random), uniform_01(random), uniform_01(random));
    }
};

}

}
