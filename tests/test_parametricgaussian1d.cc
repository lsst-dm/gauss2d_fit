#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN

#include "doctest.h"

#include <cmath>
#include <memory>

#include "parameters.h"
#include "parametricgaussian1d.h"

namespace g2f = gauss2d::fit;

TEST_CASE("ParametricGaussian1d") {
    CHECK(g2f::ParametricGaussian1D().str().size() > 0);
    auto mean_param = std::make_shared<g2f::MeanParameterD>();
    auto stddev_param = std::make_shared<g2f::StdDevParameterD>();
    double mean = mean_param->get_value();
    double stddev = stddev_param->get_value();
    auto prior = g2f::ParametricGaussian1D(mean_param, stddev_param);

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
}
