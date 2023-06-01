#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN

#include "doctest.h"

#include <cmath>
#include <memory>

#include "gaussianparametricellipse.h"
#include "parameters.h"
#include "parametricgaussian1d.h"
#include "shapeprior.h"
#include "transforms.h"

namespace g2 = gauss2d;
namespace g2f = gauss2d::fit;

TEST_CASE("ShapePrior") {
    CHECK_THROWS_AS(g2f::ShapePrior(nullptr), std::invalid_argument);
    auto ellipse = std::make_shared<g2f::GaussianParametricEllipse>();
    auto options = std::make_shared<g2f::ShapePriorOptions>();
    CHECK(options->repr().size() > 0);
    CHECK(options->repr(true).size() > 0);
    CHECK(options->str().size() > 0);
    CHECK(g2f::ShapePrior(ellipse).str().size() > 0);
    CHECK(g2f::ShapePrior(ellipse, nullptr, nullptr, options).repr().size() > 0);
    double mean_size = 2.1;
    double stddev_size = 0.7;
    double mean_axrat = 0.7;
    double stddev_axrat = 0.5;
    auto limits_axrat_logit = std::make_shared<parameters::Limits<double>>(-1e-10, 1 + 1e-10);
    auto transform_axrat = std::make_shared<g2f::LogitLimitedTransform>(limits_axrat_logit);

    auto prior_size = std::make_shared<g2f::ParametricGaussian1D>(
            std::make_shared<g2f::MeanParameter>(mean_size, nullptr,
                                                 g2f::get_transform_default<g2f::Log10Transform>()),
            std::make_shared<g2f::StdDevParameter>(stddev_size));
    auto prior_axrat = std::make_shared<g2f::ParametricGaussian1D>(
            std::make_shared<g2f::MeanParameter>(mean_axrat, nullptr, transform_axrat),
            std::make_shared<g2f::StdDevParameter>(stddev_axrat));
    auto prior = g2f::ShapePrior(ellipse, prior_size, prior_axrat, options);

    CHECK(prior.evaluate().loglike != 0);
    auto terms = prior.get_loglike_const_terms();
    CHECK(terms.size() == 2);
    CHECK(prior.repr().size() > 0);
    CHECK(prior.repr(true).size() > 0);
    CHECK(prior.str().size() > 0);
}
