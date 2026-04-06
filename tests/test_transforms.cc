#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN

#include "doctest.h"

#include <cmath>
#include <memory>

#include "lsst/gauss2d/fit/transforms.h"

namespace g2f = lsst::gauss2d::fit;

TEST_CASE("Log10Transform") {
    auto transform = g2f::Log10TransformD();
    CHECK_EQ(transform.str(), "Log10TransformD()");
    CHECK_EQ(transform.forward(10.), 1.);
    CHECK_EQ(transform.reverse(1.), 10.);
}

TEST_CASE("LogTransform") {
    auto transform = g2f::LogTransformD();
    CHECK_EQ(transform.str(), "LogTransformD()");
    CHECK_EQ(transform.forward(1.), 0.);
    CHECK_EQ(transform.reverse(0.), 1.);
}
