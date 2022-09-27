#include "sersicmix.h"
#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN

#include "doctest.h"

#include "linearsersicmixinterpolator.h"

namespace g2f = gauss2d::fit;

TEST_CASE("LinearSersicMixInterpolator")
{
    std::vector<unsigned short> orders = {4, 8};
    for(const auto order : orders) {
        auto interpolator = g2f::LinearSersicMixInterpolator(order);
        CHECK(interpolator.sersicindex_max > interpolator.sersicindex_min);
        CHECK_THROWS_AS(interpolator.get_integralsizes(0.499), std::invalid_argument);
        CHECK_THROWS_AS(interpolator.get_integralsizes(interpolator.sersicindex_max + 1e-10), std::invalid_argument);

        auto result = interpolator.get_integralsizes((interpolator.sersicindex_min + interpolator.sersicindex_max)/2.0);
        CHECK(result.size() == order);
        for(size_t ord = 0; ord < order; ++ord) {
            CHECK(result[ord].integral >= 0);
            CHECK(result[ord].sigma >= 0);
        }
        auto knots = gauss2d::fit::get_sersic_mix_knots(order);
        double sersicindices[2] = {interpolator.sersicindex_min, interpolator.sersicindex_max};
        g2f::SersicMixValues expected[2] = {knots[0], knots.back()};
        for(size_t i = 0; i < 2; ++i) {
            double sersicindex = sersicindices[i];
            const auto & knots = expected[i];
            CHECK(knots.sersicindex == sersicindex);
            auto interp = interpolator.get_integralsizes(sersicindex);
            for(size_t ord = 0; ord < order; ++ord) {
                CHECK(knots.values[ord].integral == interp[ord].integral);
                CHECK(knots.values[ord].sigma == interp[ord].sigma);
            }
        }
    }
}