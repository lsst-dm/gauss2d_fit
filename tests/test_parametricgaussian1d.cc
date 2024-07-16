#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN

#include "doctest.h"

#include <memory>

#include "lsst/gauss2d/fit/parameters.h"
#include "lsst/gauss2d/fit/parametricgaussian1d.h"

namespace g2f = lsst::gauss2d::fit;

TEST_CASE("ParametricGaussian1d") {
    CHECK_GT(g2f::ParametricGaussian1D().str().size(), 0);
    auto mean_param = std::make_shared<g2f::MeanParameterD>();
    auto stddev_param = std::make_shared<g2f::StdDevParameterD>();
    double mean = mean_param->get_value();
    double stddev = stddev_param->get_value();
    auto prior = g2f::ParametricGaussian1D(mean_param, stddev_param);

    CHECK_GT(prior.repr().size(), 0);
    CHECK_GT(prior.repr(true).size(), 0);
    CHECK_GT(prior.str().size(), 0);
    CHECK_EQ(prior.get_mean(), mean);
    mean += 1;
    prior.set_mean(mean);
    CHECK_EQ(prior.get_mean(), mean);
    CHECK_EQ(prior.get_stddev(), stddev);
    stddev += 1;
    prior.set_stddev(stddev);
    CHECK_EQ(prior.get_stddev(), stddev);
}
