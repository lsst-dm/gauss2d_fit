#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN

#include "doctest.h"

#include "lsst/gauss2d/fit/channel.h"
#include "lsst/gauss2d/fit/linearintegralmodel.h"
#include "lsst/gauss2d/fit/param_defs.h"
#include "lsst/gauss2d/fit/param_filter.h"

namespace g2f = lsst::gauss2d::fit;

TEST_CASE("LinearIntegralModel") {
    const auto c1 = g2f::Channel::make("1");
    const auto c2 = g2f::Channel::make("2");
    const auto c3 = g2f::Channel::make("3");

    auto p1 = std::make_shared<g2f::IntegralParameterD>(1);
    auto p2 = std::make_shared<g2f::IntegralParameterD>(2);

    g2f::LinearIntegralModel::Data data = {{*c1, p1}, {*c2, p2}};

    auto integral = std::make_shared<g2f::LinearIntegralModel>(&data);
    CHECK_GT(integral->str().size(), 0);
    CHECK_EQ(integral->size(), 2);

    const auto channels = integral->get_channels();
    CHECK_EQ(*channels.begin(), *c1);

    integral->at(*c1)->set_value(0);

    g2f::ParamCRefs params;
    integral->get_parameters_const(params);
    CHECK_EQ(params.size(), 2);
    CHECK_EQ(params.at(0), *p1);

    g2f::ParamFilter filter{true, true, true, true};
    g2f::ParamRefs params2;
    integral->get_parameters(params2, &filter);
    CHECK_EQ(params2.size(), 2);
    filter.channel = *c2;
    integral->get_parameters(params2, &filter);
    CHECK_EQ(params2.size(), 3);
    CHECK_EQ(params2.at(2), *p2);
    filter.channel = *c3;
    integral->get_parameters(params2, &filter);
    CHECK_EQ(params2.size(), 3);
}

TEST_CASE("LinearIntegralModel(default)") {
    const auto& channel = g2f::Channel::NONE();
    auto integral = std::make_shared<g2f::LinearIntegralModel>(nullptr);
    double value = 2.;
    integral->at(channel)->set_value(value);
    CHECK_EQ(integral->get_integral(channel), value);
}