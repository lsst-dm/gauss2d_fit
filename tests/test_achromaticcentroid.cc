#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN

#include "doctest.h"

#include "lsst/gauss2d/fit/achromaticcentroid.h"

namespace g2f = lsst::gauss2d::fit;

TEST_CASE("ChromaticCentroid") {
    const unsigned int DIM = 20;

    g2f::AchromaticCentroid cen_default{};
    CHECK_EQ(cen_default.n_channels(), 0);
    CHECK_EQ(cen_default.get_channels().size(), 0);

    auto cenx = std::make_shared<g2f::CentroidXParameterD>(DIM / 2.);
    auto ceny = std::make_shared<g2f::CentroidYParameterD>(DIM / 2.);

    auto cens = std::make_shared<g2f::CentroidParameters>(cenx, ceny);

    auto centroid = std::make_shared<g2f::AchromaticCentroid>(cens);
    auto cen_at = centroid->at(g2f::Channel::NONE());

    CHECK_EQ(cen_at, cens);

    g2f::ParamCRefs params{};
    centroid->get_parameters_const(params);
    CHECK_EQ(params.size(), 2);
    auto params_uniq = g2f::nonconsecutive_unique(params);
    CHECK_EQ(params_uniq.size(), 2);
}