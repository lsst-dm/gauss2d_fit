#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN

#include "doctest.h"

#include "lsst/gauss2d/fit/channel.h"
#include "lsst/gauss2d/fit/fractionalintegralmodel.h"
#include "lsst/gauss2d/fit/gaussiancomponent.h"
#include "lsst/gauss2d/fit/linearintegralmodel.h"
#include "lsst/gauss2d/fit/param_defs.h"
#include "lsst/gauss2d/fit/source.h"

namespace g2f = lsst::gauss2d::fit;

TEST_CASE("FractionalIntegralModel") {
    const auto c1 = g2f::Channel::make("1");
    const auto c2 = g2f::Channel::make("2");
    const auto c3 = g2f::Channel::make("3");

    auto p1 = std::make_shared<g2f::IntegralParameterD>(1);
    auto p2 = std::make_shared<g2f::IntegralParameterD>(2);
    auto p3 = std::make_shared<g2f::IntegralParameterD>(3);

    g2f::LinearIntegralModel::Data data = {{*c1, p1}, {*c2, p2}, {*c3, p3}};
    auto integral = std::make_shared<g2f::LinearIntegralModel>(&data);

    const auto& channel = *c1;

    g2f::FractionalIntegralModel::Data data_frac{{*c1, std::make_shared<g2f::ProperFractionParameterD>(0.5)},
                                                 {*c2, std::make_shared<g2f::ProperFractionParameterD>(0.25)},
                                                 {*c3, std::make_shared<g2f::ProperFractionParameterD>(0.1)}};
    auto frac = g2f::FractionalIntegralModel::make(data_frac, integral);

    // Check that the parent is set correctly and gives the right value
    CHECK_EQ(&(frac->get_parent_model()), integral.get());
    CHECK_EQ(frac->at(channel)->get_value(), 0.5);
    CHECK_EQ(frac->get_integral_remainder(channel), 0.5);

    // Set the fraction to 100% (previously 0)
    frac->at(channel)->set_value(1.);
    CHECK_EQ(frac->get_integral(channel), 1);
    CHECK_EQ(frac->get_integral_remainder(channel), 0);

    // Check that setting to 25% works
    double frac1 = 0.25;
    frac->at(channel)->set_value(frac1);
    CHECK_EQ(frac->at(channel)->get_value(), frac1);
    CHECK_EQ(frac->get_integral(channel), frac1);
    CHECK_EQ(frac->get_integral_remainder(channel), (1. - frac1));

    // Check that the frac gives both parameters (parent's integral and its own fraction)
    g2f::ParamCRefs params{};
    CHECK_EQ(frac->get_parameters_const(params).size(), 6);

    // Add another frac and validate it
    g2f::FractionalIntegralModel::Data data_frac2{
            {*c1, std::make_shared<g2f::ProperFractionParameterD>(0.5)},
            {*c2, std::make_shared<g2f::ProperFractionParameterD>(0.5)},
            {*c3, std::make_shared<g2f::ProperFractionParameterD>(0.5)}};
    auto frac2 = g2f::FractionalIntegralModel::make(data_frac2, frac);
    frac2->at(channel)->set_value(1.);
    CHECK_EQ(frac2->at(channel)->get_value(), 1);
    CHECK_EQ(frac2->get_integral(channel), (1 - frac1));
    CHECK_EQ(frac->find_model(*integral), nullptr);

    std::shared_ptr<g2f::IntegralModel> integral_base = std::make_shared<g2f::LinearIntegralModel>(nullptr);
    CHECK_EQ(std::dynamic_pointer_cast<g2f::FractionalIntegralModel>(integral_base), nullptr);

    g2f::FractionalIntegralModel::Data data_frac3{
            {*c1, std::make_shared<g2f::ProperFractionParameterD>(1.0, nullptr, nullptr, nullptr, true)},
            {*c2, std::make_shared<g2f::ProperFractionParameterD>(1.0, nullptr, nullptr, nullptr, true)},
            {*c3, std::make_shared<g2f::ProperFractionParameterD>(1.0)}};
    // All prop fraction parameters must be fixed
    CHECK_THROWS_AS(g2f::FractionalIntegralModel::make(data_frac3, frac2, true), std::invalid_argument);
    data_frac3[1].second->set_value(0.0);
    data_frac3[2].second->set_fixed(true);
    // Still no good, frac must be 1.0
    CHECK_THROWS_AS(g2f::FractionalIntegralModel::make(data_frac3, frac2, true), std::invalid_argument);
    data_frac3[1].second->set_value(1.0);
    auto last = g2f::FractionalIntegralModel::make(data_frac3, frac2, true);

    // Check that interactions with GaussianComponents work
    auto comp1 = std::make_shared<g2f::GaussianComponent>(nullptr, nullptr, frac);
    auto comp2 = std::make_shared<g2f::GaussianComponent>(nullptr, nullptr, frac2);
    auto comp3 = std::make_shared<g2f::GaussianComponent>(nullptr, nullptr, last);
    g2f::Components comps{comp1, comp2, comp3};

    auto source = std::make_shared<g2f::Source>(comps);
    auto gaussians = source->get_gaussians(*c1);
    auto g1 = gaussians->at(0);
    CHECK_EQ(g1.get_integral_value(), frac->at(channel)->get_value());
    // TODO: Add more checks
}
