#include "sersicmixcomponent.h"
#include "parameters.h"
#include <stdexcept>
#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN

#include "doctest.h"

#include "channel.h"
#include "linearintegralmodel.h"
#include "sersicmixcomponent.h"

namespace g2f = gauss2d::fit;

TEST_CASE("SersicMixComponentIndexParameter")
{
    auto param = g2f::SersicMixComponentIndexParameter(0.5);
    CHECK_THROWS(param.set_limits(std::make_shared<const parameters::Limits<double>>(0.3, 8.1)));
    CHECK(param.get_limits_maximal().get_min() == 0.5);
    CHECK(param.get_limits().get_min() == 0.5);
    CHECK_THROWS_AS(param.set_value(0.5-1e-10), std::runtime_error);
    param.set_limits(std::make_shared<const parameters::Limits<double>>(1.0, 4.0));
    CHECK_THROWS_AS(param.set_value(1.0 - 1e-10), std::runtime_error);
    CHECK_THROWS_AS(param.set_value(4.0 + 1e-10), std::runtime_error);
}

TEST_CASE("SersicMixComponent")
{
    auto sersic_n = std::make_shared<g2f::SersicMixComponentIndexParameter>(1.0);
    CHECK(sersic_n->order == 4);
    
    const auto & C = g2f::Channel::NONE();
    
    auto reff_x = std::make_shared<g2f::ReffXParameter>(1.0);
    auto reff_y = std::make_shared<g2f::ReffYParameter>(1.0);

    const double total = 2.;

    g2f::LinearIntegralModel::Data data = {{C, std::make_shared<g2f::IntegralParameter>(total)}};

    auto comp = std::make_shared<g2f::SersicMixComponent>(
        std::make_shared<g2f::SersicParametricEllipse>(reff_x, reff_y),
        std::make_shared<g2f::CentroidParameters>(),
        std::make_shared<g2f::LinearIntegralModel>(&data),
        sersic_n
    );
    CHECK(comp->str().size() > 0);

    CHECK(comp->get_n_gaussians(C) == sersic_n->order);
    const auto gaussians = comp->get_gaussians(C);
    CHECK(gaussians->size() == sersic_n->order);
    
    double integral = 0;
    double sigma_old = 1e31;
    for(const auto & g : *gaussians) {
        integral += g->get_integral_value();
        const auto & ell = g->get_ellipse_const();
        double sigma = ell.get_sigma_x();
        CHECK(sigma > 0);
        CHECK(sigma < sigma_old);
        CHECK(ell.get_sigma_y() == sigma);
        CHECK(ell.get_rho() == 0);
        sigma_old = sigma;
    }
    CHECK(integral == total);

    sersic_n->set_value(0.5);
    const auto gauss = comp->get_gaussians(C);
    CHECK(gauss->at(0).get_integral_value() == total);
    for(size_t i = 1; i < sersic_n->order; ++i) CHECK(gauss->at(i).get_integral_value() == 0);

    // 3 ellipse, 2 centroid, one integral, one Sersic index 
    const size_t n_params = 7;

    g2f::ParamRefs params;
    CHECK(comp->get_parameters(params).size() == n_params);

    g2f::ParamCRefs cparams;
    CHECK(comp->get_parameters_const(cparams).size() == n_params);
}
