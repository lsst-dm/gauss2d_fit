#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN

#include "doctest.h"

#include "channel.h"
#include "linearintegralmodel.h"
#include "param_defs.h"

namespace g2f = gauss2d::fit;

TEST_CASE("LinearIntegralModel")
{
    const auto c1 = g2f::Channel::make("1");
    const auto c2 = g2f::Channel::make("2");

    auto p1 = std::make_shared<g2f::IntegralParameter>(1);
    auto p2 = std::make_shared<g2f::IntegralParameter>(2);

    g2f::LinearIntegralModel::Data data = {{*c1, p1}, {*c2, p2}};

    auto integral = std::make_shared<g2f::LinearIntegralModel>(&data);
    CHECK(integral->str().size() > 0);
    CHECK(integral->size() == 2);

    const auto channels = integral->get_channels();
    CHECK(*channels.begin() == *c1);

    integral->at(*c1)->set_value(0);

    g2f::ParamCRefs params;
    integral->get_parameters_const(params);
    CHECK(params.size() == 2);
    CHECK(params.at(0) == *p1);
}