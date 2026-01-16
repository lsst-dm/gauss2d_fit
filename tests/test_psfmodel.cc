#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN

#include "doctest.h"

#include <memory>

#include "lsst/gauss2d/fit/channel.h"
#include "lsst/gauss2d/fit/component.h"
#include "lsst/gauss2d/fit/gaussiancomponent.h"
#include "lsst/gauss2d/fit/psfmodel.h"

namespace g2f = lsst::gauss2d::fit;

const auto& CHANNEL_NONE = g2f::Channel::NONE();

TEST_CASE("PsfModel") {
    auto comp = std::make_shared<g2f::GaussianComponent>();
    const auto& comp_const = *comp;
    std::vector<std::shared_ptr<g2f::Component>> comps;
    comps.push_back(std::move(comp));
    comps.push_back(std::make_shared<g2f::GaussianComponent>());

    auto psfmodel = std::make_shared<g2f::PsfModel>(comps);

    g2f::ParamCRefs params{};
    CHECK_EQ(psfmodel->get_components().size(), 2);
    // 2 comps x (2 centroid, 3 ellipse, 1 integral)
    CHECK_EQ(psfmodel->get_parameters_const(params).size(), 12);
    CHECK_EQ(psfmodel->get_n_gaussians(CHANNEL_NONE), 2);
    auto gaussians = psfmodel->get_gaussians(CHANNEL_NONE);
    CHECK_EQ(gaussians->size(), 2);
    const auto& g0 = gaussians->at(0);
    const auto gaussians_comp = comp_const.get_gaussians(CHANNEL_NONE);
    CHECK_EQ(gaussians_comp->size(), 1);
    const auto& c0 = gaussians_comp->at(0);
    CHECK_EQ(g0.get_centroid_const(), c0.get_centroid_const());
    CHECK_EQ(g0.get_ellipse_const(), c0.get_ellipse_const());
    CHECK_EQ(g0.get_integral_const(), c0.get_integral_const());
    CHECK_EQ(g0, c0);
}

TEST_CASE("PsfModelNontrivialComponent") {
    g2f::LinearIntegralModel::Data data = {{CHANNEL_NONE, std::make_shared<g2f::IntegralParameterD>(1)}};

    auto comp = std::make_shared<g2f::GaussianComponent>(
        std::make_shared<g2f::GaussianParametricEllipse>(1.0, 1.0, 0.0),
        nullptr,
        std::make_shared<g2f::LinearIntegralModel>(&data)
    );
    std::vector<std::shared_ptr<g2f::Component>> comps;
    comps.push_back(std::move(comp));
    comps.push_back(std::make_shared<g2f::GaussianComponent>());

    auto psfmodel = std::make_shared<g2f::PsfModel>(comps);
    CHECK_EQ(psfmodel->get_components().size(), 2);
}

