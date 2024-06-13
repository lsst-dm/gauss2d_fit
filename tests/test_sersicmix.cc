#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN

#include "doctest.h"

#include "lsst/gauss2d/fit/sersicmix.h"

namespace g2f = lsst::gauss2d::fit;

TEST_CASE("SersicMixtures") {
    std::vector<unsigned short> orders = {4, 8};
    for (const auto order : orders) {
        const auto knots = g2f::get_sersic_mix_knots_copy(order);
        CHECK_EQ(knots[0].sersicindex, 0.5);
        CHECK_GT(knots.back().sersicindex, 1);
    }
}