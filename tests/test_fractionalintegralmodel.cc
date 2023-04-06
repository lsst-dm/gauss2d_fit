#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN

#include "doctest.h"

#include "channel.h"
#include "fractionalintegralmodel.h"
#include "linearintegralmodel.h"
#include "param_defs.h"

namespace g2f = gauss2d::fit;

TEST_CASE("FractionalIntegralModel") {
    const auto& channel = g2f::Channel::NONE();
    auto integral = std::make_shared<g2f::LinearIntegralModel>(nullptr);
    double value = 2.;
    integral->at(channel)->set_value(value);
    auto frac = g2f::FractionalIntegralModel::make(std::nullopt, integral);

    // Check that the parent is set correctly and gives the right value
    CHECK(&(frac->get_parent_model()) == integral.get());
    CHECK(frac->at(channel)->get_value() == 0);
    CHECK(frac->get_integral_remainder(channel) == value);

    // Set the fraction to 100% (previously 0)
    frac->at(channel)->set_value(1.);
    CHECK(frac->get_integral(channel) == value);
    CHECK(frac->get_integral_remainder(channel) == 0);

    // Check that setting to 25% works
    double frac1 = 0.25;
    frac->at(channel)->set_value(frac1);
    CHECK(frac->at(channel)->get_value() == frac1);
    CHECK(frac->get_integral(channel) == frac1 * value);
    CHECK(frac->get_integral_remainder(channel) == (1. - frac1) * value);

    // Check that the frac gives both parameters (parent's integral and its own fraction)
    g2f::ParamCRefs params{};
    CHECK(frac->get_parameters_const(params).size() == 2);

    // Add another frac and valifate it
    auto frac2 = g2f::FractionalIntegralModel::make(std::nullopt, frac);
    frac2->at(channel)->set_value(1.);
    CHECK(frac2->at(channel)->get_value() == 1.);
    CHECK(frac2->get_integral(channel) == (1 - frac1) * value);
    CHECK(frac->find_model(*integral) == nullptr);

    std::shared_ptr<g2f::IntegralModel> integral_base = std::make_shared<g2f::LinearIntegralModel>(nullptr);
    CHECK(std::dynamic_pointer_cast<g2f::FractionalIntegralModel>(integral_base) == nullptr);
}