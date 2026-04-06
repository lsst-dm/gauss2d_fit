#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN

#include "doctest.h"

#include "lsst/gauss2d/fit/parameters.h"

namespace g2f = lsst::gauss2d::fit;

TEST_CASE("SersicParameter") {
    auto param = g2f::SersicIndexParameterD();
    CHECK_EQ(param.str(), "SersicIndexParameterD(value=" + std::to_string(param.get_value()) + ", " + ")");
}
