#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN

#include "doctest.h"

#include "lsst/gauss2d/ellipse.h"

#include "lsst/gauss2d/fit/sersicmix.h"

namespace g2d = lsst::gauss2d;
namespace g2f = lsst::gauss2d::fit;

TEST_CASE("SersicMixtures") {
    std::vector<unsigned short> orders = {4, 8};
    for (const auto order : orders) {
        const auto knots = g2f::get_sersic_mix_knots_copy(order);
        double sersic_index = 0;
        // Can't make n<0.5 out of concentric Gaussians, so 0.5 should be first
        CHECK_EQ(knots[0].sersicindex, 0.5);
        const auto& knot_values_0 = knots[0].values;
        CHECK_EQ(knot_values_0[0].integral, 1);
        CHECK_EQ(knot_values_0[0].sigma, lsst::gauss2d::M_SIGMA_HWHM);
        for (size_t i = 1; i < knot_values_0.size(); ++i) {
            CHECK_EQ(knot_values_0[i].integral, 0);
        }
        // Check that it goes beyond n=4 at least
        CHECK_GT(knots.back().sersicindex, 4);

        for (const auto& knot : knots) {
            // Check that the Sersic index increases
            CHECK_GT(knot.sersicindex, sersic_index);
            sersic_index = knot.sersicindex;
            double integral = 0;
            double radius = std::numeric_limits<double>::infinity();
            for (const auto& knot_values : knot.values) {
                integral += knot_values.integral;
                // Check that the radius decreases
                CHECK_LT(knot_values.sigma, radius);
                radius = knot_values.sigma;
            }
            // Check for normalization
            CHECK_LT(std::abs(1 - integral), 1e-14);
        }
    }
}