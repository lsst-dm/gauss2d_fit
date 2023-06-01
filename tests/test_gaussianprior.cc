#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN

#include "doctest.h"

#include <cmath>
#include <memory>

#include "gaussianprior.h"
#include "parameters.h"

namespace g2 = gauss2d;
namespace g2f = gauss2d::fit;

TEST_CASE("GaussianPrior") {
    auto param = std::make_shared<g2f::IntegralParameter>();
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
    prior.set_transformed(true);
    CHECK(prior.get_transformed() == true);
}
