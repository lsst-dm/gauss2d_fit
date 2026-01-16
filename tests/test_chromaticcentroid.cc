#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN

#include "doctest.h"

#include "lsst/gauss2d/fit/chromaticcentroid.h"
#include "lsst/gauss2d/fit/util.h"

namespace g2f = lsst::gauss2d::fit;

TEST_CASE("ChromaticCentroid") {
    const unsigned int DIM = 20;

    g2f::ChromaticCentroid cen_default{};
    CHECK_EQ(cen_default.n_channels(), 1);
    CHECK_EQ(cen_default.get_channels()[0], g2f::Channel::NONE());

    auto cenx = std::make_shared<g2f::CentroidXParameterD>(DIM / 2.);
    auto ceny = std::make_shared<g2f::CentroidYParameterD>(DIM / 2.);

    auto cen1 = std::make_shared<g2f::CentroidParameters>(cenx, ceny);
    auto cen2 = std::make_shared<g2f::CentroidParameters>(cenx, ceny);

    auto cenx3 = std::make_shared<g2f::CentroidXParameterD>(DIM / 2. + 0.1);
    auto ceny3 = std::make_shared<g2f::CentroidYParameterD>(DIM / 2. + 0.1);

    auto cen3 = std::make_shared<g2f::CentroidParameters>(cenx3, ceny3);

    const auto c1 = g2f::Channel::make("1");
    const auto c2 = g2f::Channel::make("2");
    const auto c3 = g2f::Channel::make("3");
    const auto c4 = g2f::Channel::make("4");

    g2f::ChromaticCentroid::Data data = {{*c1, cen1}, {*c2, cen2}, {*c3, cen3}, {*c4, cen1}};

    auto centroid = std::make_shared<g2f::ChromaticCentroid>(data);
    auto cen_c1 = centroid->at(*c1);

    CHECK_EQ(cen_c1, cen1);
    CHECK_EQ(centroid->at(*c4), cen1);
    CHECK_NE(centroid->at(*c3), cen1);
    CHECK_NE(centroid->at(*c2), centroid->at(*c3));

    cenx->set_value(0.);
    CHECK_EQ(centroid->at(*c4)->get_x(), 0);

    g2f::ParamCRefs params{};
    centroid->get_parameters_const(params);
    CHECK_EQ(params.size(), 8);
    auto params_uniq = g2f::nonconsecutive_unique(params);
    CHECK_EQ(params_uniq.size(), 4);
}