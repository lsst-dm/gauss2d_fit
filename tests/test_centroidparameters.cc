#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN

#include "doctest.h"

#include "centroidparameters.h"

namespace g2f = gauss2d::fit;

TEST_CASE("CentroidParameters")
{
    const unsigned int DIM = 20;

    auto cenx = std::make_shared<g2f::CentroidXParameter>(DIM/2.);
    auto ceny = std::make_shared<g2f::CentroidYParameter>(DIM/2.);

    auto cen = std::make_shared<g2f::CentroidParameters>(cenx, ceny);
    CHECK(cen->get_xy() == std::array<double, 2>{DIM/2., DIM/2.});

    g2f::ParamCRefs params{};
    CHECK(cen->get_parameters_const(params).size() == 2);
}