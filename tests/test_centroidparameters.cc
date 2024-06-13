#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN

#include "doctest.h"

#include "lsst/gauss2d/fit/centroidparameters.h"

namespace g2f = lsst::gauss2d::fit;

TEST_CASE("CentroidParameters") {
    const unsigned int DIM = 20;

    auto cenx = std::make_shared<g2f::CentroidXParameterD>(DIM / 2.);
    auto ceny = std::make_shared<g2f::CentroidYParameterD>(DIM / 2.);

    auto cen = std::make_shared<g2f::CentroidParameters>(cenx, ceny);
    CHECK_EQ(cen->get_xy(), std::array<double, 2>{DIM / 2., DIM / 2.});

    g2f::ParamCRefs params{};
    CHECK_EQ(cen->get_parameters_const(params).size(), 2);
}