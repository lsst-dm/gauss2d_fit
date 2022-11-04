#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN

#include "doctest.h"

#include <memory>

#include "channel.h"
#include "component.h"
#include "gaussiancomponent.h"
#include "psfmodel.h"

namespace g2f = gauss2d::fit;

TEST_CASE("PsfModel") {
    const auto & C = g2f::Channel::NONE();
    auto comp = std::make_shared<g2f::GaussianComponent>();
    const auto & comp_const = *comp;
    std::vector<std::shared_ptr<g2f::Component>> comps;
    comps.push_back(std::move(comp));
    comps.push_back(std::make_shared<g2f::GaussianComponent>());

    auto psfmodel = std::make_shared<g2f::PsfModel>(comps);

    g2f::ParamCRefs params{};
    CHECK(psfmodel->get_components().size() == 2);
    // 2 comps x (2 centroid, 3 ellipse, 1 integral)
    CHECK(psfmodel->get_parameters_const(params).size() == 12);
    CHECK(psfmodel->get_n_gaussians(C) == 2);
    auto gaussians = psfmodel->get_gaussians(C);
    CHECK(gaussians->size() == 2);
    const auto & g0 = gaussians->at(0);
    const auto gaussians_comp = comp_const.get_gaussians(C);
    CHECK(gaussians_comp->size() == 1);
    const auto & c0 = gaussians_comp->at(0);
    CHECK(g0.get_centroid_const() == c0.get_centroid_const());
    CHECK(g0.get_ellipse_const() == c0.get_ellipse_const());
    CHECK(g0.get_integral_const() == c0.get_integral_const());
    CHECK(g0 == c0);
}
