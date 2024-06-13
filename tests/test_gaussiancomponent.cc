#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN

#include "doctest.h"

#include "lsst/gauss2d/fit/channel.h"
#include "lsst/gauss2d/fit/gaussiancomponent.h"

namespace g2d = lsst::gauss2d;
namespace g2f = lsst::gauss2d::fit;

TEST_CASE("GaussianComponent") {
    const auto& C = g2f::Channel::NONE();
    auto comp = std::make_shared<g2f::GaussianComponent>();
    CHECK_GT(comp->str().size(), 0);

    g2f::ParamCRefs params{};
    // 2 centroid, 3 ellipse, 1 integral
    CHECK_EQ(comp->get_parameters_const(params).size(), 6);

    CHECK_EQ(comp->get_n_gaussians(C), 1);
    auto gaussians = comp->get_gaussians(C);
    CHECK_EQ(gaussians->size(), 1);
    const auto& g0 = gaussians->at(0);
    const g2d::Gaussian g1{};
    CHECK_EQ(g0.get_centroid_const(), g1.get_centroid_const());
    CHECK_EQ(g0.get_ellipse_const(), g1.get_ellipse_const());
    CHECK_EQ(g0.get_integral_const(), g1.get_integral_const());

    CHECK_EQ(gaussians->at_const(0), g1);

    auto gaussians_default = g2f::GaussianComponent::make_uniq_default_gaussians({1., 3.});
    CHECK_EQ(gaussians_default.size(), 2);
}
