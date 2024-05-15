#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN

#include "doctest.h"

#include <cmath>
#include <memory>

#include "gaussianprior.h"
#include "parameters.h"
#include "transforms.h"
#include "util.h"

namespace g2 = gauss2d;
namespace g2f = gauss2d::fit;

TEST_CASE("GaussianPrior") {
    auto param = std::make_shared<g2f::IntegralParameterD>();
    double mean = param->get_value();
    double stddev = 1;
    auto prior = g2f::GaussianPrior(param, mean, stddev, false);

    CHECK(prior.evaluate().loglike == 0);
    auto terms = prior.get_loglike_const_terms();
    CHECK(terms.size() == 1);
    CHECK(terms[0] == -log(sqrt(2 * M_PI)));
    CHECK(prior.repr().size() > 0);
    CHECK(prior.repr(true).size() > 0);
    CHECK(prior.str().size() > 0);
    CHECK(prior.get_mean() == mean);
    mean += 1;
    prior.set_mean(mean);
    CHECK(prior.get_mean() == mean);
    CHECK(prior.get_stddev() == stddev);
    stddev += 1;
    prior.set_stddev(stddev);
    CHECK(prior.get_stddev() == stddev);
    CHECK(prior.get_transformed() == false);

    param = std::make_shared<g2f::IntegralParameterD>(1.0, nullptr, std::make_shared<g2f::Log10Transform>());
    prior = g2f::GaussianPrior(param, param->get_value_transformed(), stddev, true);
    CHECK(prior.get_transformed() == true);
    auto eval = prior.evaluate(true);
    CHECK(eval.residuals[0] == 0);
    CHECK(eval.jacobians.at(*param)[0] == 1.0 / stddev);

    const double delta = 1e-6;
    std::vector<double> test_x{-0.34, 0., 1.21};
    for (double x : test_x) {
        param->set_value_transformed(x - delta / 2.);
        auto eval = prior.evaluate(true);
        param->set_value_transformed(x + delta / 2.);
        const double dll_dx = (prior.evaluate().loglike - eval.loglike) / delta;
        param->set_value_transformed(x);
        CHECK(g2f::isclose(dll_dx, eval.compute_dloglike_dx(*param), delta, delta).isclose);
    }
}
