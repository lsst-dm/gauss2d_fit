#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN

#include "doctest.h"

#include <memory>

#include "lsst/gauss2d/fit/channel.h"
#include "lsst/gauss2d/fit/component.h"
#include "lsst/gauss2d/fit/gaussiancomponent.h"
#include "lsst/gauss2d/fit/source.h"

namespace g2f = lsst::gauss2d::fit;

TEST_CASE("Source") {
    const auto& C = g2f::Channel::NONE();
    auto comp = std::make_unique<g2f::GaussianComponent>();
    const auto& comp_const = *comp;
    std::vector<std::shared_ptr<g2f::Component>> comps;
    comps.push_back(std::move(comp));
    comps.push_back(std::make_unique<g2f::GaussianComponent>());
    auto source = std::make_shared<g2f::Source>(comps);

    g2f::ParamCRefs params{};
    CHECK_EQ(source->get_components().size(), 2);
    // 2 comps x (2 centroid, 3 ellipse, 1 integral)
    CHECK_EQ(source->get_parameters_const(params).size(), 12);

    CHECK_EQ(source->get_n_gaussians(C), 2);
    auto gaussians = source->get_gaussians(C);
    CHECK_EQ(gaussians->size(), 2);
    const auto& g0 = gaussians->at(0);
    const auto& gaussians_comp = comp_const.get_gaussians(C);
    CHECK_EQ(gaussians_comp->size(), 1);
    const auto& c0 = gaussians_comp->at(0);
    CHECK_EQ(g0.get_centroid_const(), c0.get_centroid_const());
    CHECK_EQ(g0.get_ellipse_const(), c0.get_ellipse_const());
    CHECK_EQ(g0.get_integral_const(), c0.get_integral_const());
    CHECK_EQ(g0, c0);
}
