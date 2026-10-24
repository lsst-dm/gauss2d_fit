#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN

#include "doctest.h"

#include "lsst/gauss2d/fit/gaussianparametricellipse.h"

namespace g2d = lsst::gauss2d;
namespace g2f = lsst::gauss2d::fit;

TEST_CASE("GaussianParametricEllipse") {
    const double SIZE = 5.;
    const double CORR = 0.5;

    auto sig_x = std::make_shared<g2f::SigmaXParameterD>(SIZE);
    auto sig_y = std::make_shared<g2f::SigmaYParameterD>(SIZE);
    auto rho = std::make_shared<g2f::RhoParameterD>(CORR);
    CHECK_EQ(sig_x->get_value(), SIZE);
    CHECK_EQ(sig_y->get_value(), SIZE);
    CHECK_EQ(rho->get_value(), CORR);

    auto ellipse = std::make_shared<g2f::GaussianParametricEllipse>(sig_x, sig_y, rho);

    CHECK_EQ(ellipse->get_sigma_x(), SIZE);
    CHECK_EQ(ellipse->get_sigma_y(), SIZE);
    CHECK_EQ(ellipse->get_rho(), CORR);

    CHECK_EQ(ellipse->get_hwhm_x(), SIZE * g2d::M_HWHM_SIGMA);
    CHECK_EQ(ellipse->get_hwhm_y(), SIZE * g2d::M_HWHM_SIGMA);

    ellipse->set_hwhm_x(SIZE);
    ellipse->set_hwhm_y(2 * SIZE);
    ellipse->set_rho(-CORR);

    CHECK_EQ(ellipse->get_hwhm_x(), SIZE);
    CHECK_EQ(ellipse->get_hwhm_y(), 2 * SIZE);
    CHECK_EQ(ellipse->get_sigma_x(), SIZE * g2d::M_SIGMA_HWHM);
    CHECK_EQ(ellipse->get_sigma_x(), ellipse->get_size_x());
    CHECK_EQ(ellipse->get_sigma_y(), 2 * SIZE * g2d::M_SIGMA_HWHM);
    CHECK_EQ(ellipse->get_sigma_y(), ellipse->get_size_y());
    CHECK_EQ(ellipse->get_rho(), -CORR);

    CHECK_EQ(sig_x->get_value(), SIZE * g2d::M_SIGMA_HWHM);
    CHECK_EQ(sig_y->get_value(), 2 * SIZE * g2d::M_SIGMA_HWHM);
    CHECK_EQ(rho->get_value(), -CORR);
}
