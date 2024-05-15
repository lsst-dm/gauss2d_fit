#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN

#include "doctest.h"

#include "gaussianparametricellipse.h"
// #include "parameters.h"

namespace g2f = gauss2d::fit;

TEST_CASE("GaussianParametricEllipse") {
    const double SIZE = 5.;

    auto sig_x = std::make_shared<g2f::SigmaXParameterD>(SIZE);
    auto sig_y = std::make_shared<g2f::SigmaYParameterD>(SIZE);
    auto rho = std::make_shared<g2f::RhoParameterD>(0.5);
    CHECK(sig_x->get_value() == SIZE);
    CHECK(sig_y->get_value() == SIZE);
    CHECK(rho->get_value() == 0.5);

    auto ellipse = std::make_shared<gauss2d::fit::GaussianParametricEllipse>(sig_x, sig_y, rho);
}
