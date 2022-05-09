#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN

#include "doctest.h"

#include "ellipseparameters.h"
//#include "parameters.h"

namespace g2f = gauss2d::fit;

TEST_CASE("EllipseParameters")
{
    const double SIZE = 5.;

    auto sig_x = std::make_shared<g2f::SigmaXParameter>(SIZE);
    auto sig_y = std::make_shared<g2f::SigmaYParameter>(SIZE);
    auto rho = std::make_shared<g2f::RhoParameter>(0.5);
    CHECK(sig_x->get_value() == SIZE);
    CHECK(sig_y->get_value() == SIZE);
    CHECK(rho->get_value() == 0.5);
}