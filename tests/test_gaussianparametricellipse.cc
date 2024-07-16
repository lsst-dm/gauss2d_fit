#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN

#include "doctest.h"

#include "lsst/gauss2d/fit/gaussianparametricellipse.h"

namespace g2f = lsst::gauss2d::fit;

TEST_CASE("GaussianParametricEllipse") {
    const double SIZE = 5.;

    auto sig_x = std::make_shared<g2f::SigmaXParameterD>(SIZE);
    auto sig_y = std::make_shared<g2f::SigmaYParameterD>(SIZE);
    auto rho = std::make_shared<g2f::RhoParameterD>(0.5);
    CHECK_EQ(sig_x->get_value(), SIZE);
    CHECK_EQ(sig_y->get_value(), SIZE);
    CHECK_EQ(rho->get_value(), 0.5);

    auto ellipse = std::make_shared<g2f::GaussianParametricEllipse>(sig_x, sig_y, rho);
}
