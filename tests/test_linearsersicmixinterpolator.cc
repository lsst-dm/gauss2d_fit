#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN

#include "doctest.h"

#include "lsst/gauss2d/fit/linearsersicmixinterpolator.h"
#include "lsst/gauss2d/fit/sersicmix.h"
#include "lsst/gauss2d/fit/util.h"

namespace g2f = lsst::gauss2d::fit;

TEST_CASE("LinearSersicMixInterpolator") {
    std::vector<unsigned short> orders = {4, 8};
    for (const auto order : orders) {
        auto interpolator = g2f::LinearSersicMixInterpolator(order);
        CHECK(interpolator.get_interptype() == g2f::InterpType::linear);
        CHECK(interpolator.get_sersicindex_max() > interpolator.get_sersicindex_min());
        CHECK_THROWS_AS(interpolator.get_integralsizes(0.499), std::invalid_argument);
        CHECK_THROWS_AS(interpolator.get_integralsizes(interpolator.get_sersicindex_max() + 1e-10),
                        std::invalid_argument);

        double diff = interpolator.get_sersicindex_max() - interpolator.get_sersicindex_min();

        for (const double factor : {0., 0.5, 1.}) {
            double delta = 1e-6;
            const double sersicindex = interpolator.get_sersicindex_min() + diff * factor;
            auto result = interpolator.get_integralsizes(sersicindex);
            CHECK(result.size() == order);
            for (size_t ord = 0; ord < order; ++ord) {
                CHECK(result[ord].integral >= 0);
                CHECK(result[ord].sigma >= 0);
            }

            double factor_new = factor + delta;
            if (!(factor_new <= 1)) {
                delta = -delta;
                factor_new = factor + delta;
            }
            const double sersicindex_new = interpolator.get_sersicindex_min() + diff * factor_new;
            const double dsersicindex = sersicindex_new - sersicindex;
            auto result_new = interpolator.get_integralsizes(sersicindex_new);

            auto result_derivs = interpolator.get_integralsizes_derivs(sersicindex);
            for (size_t ord = 0; ord < order; ++ord) {
                CHECK(result_new[ord].integral >= 0);
                CHECK(result_new[ord].sigma >= 0);
                double deriv = result_derivs[ord].integral;
                double findiff = (result_new[ord].integral - result[ord].integral) / dsersicindex;
                auto close_integral = g2f::isclose(deriv, findiff);
                std::string msg = "sersicindex= " + std::to_string(sersicindex) + " integral["
                                  + std::to_string(ord) + "/" + std::to_string(order)
                                  + "] deriv,findiff=" + std::to_string(deriv) + "," + std::to_string(findiff)
                                  + " (from " + std::to_string(result_new[ord].integral) + "-"
                                  + std::to_string(result[ord].integral) + " result=" + close_integral.str();
                CHECK_MESSAGE(close_integral.isclose, msg, delta, delta);
                deriv = result_derivs[ord].sigma;
                findiff = (result_new[ord].sigma - result[ord].sigma) / dsersicindex;
                auto close_sigma = g2f::isclose(deriv, findiff);
                msg = "sigma[" + std::to_string(ord) + "] deriv,findiff=" + std::to_string(deriv) + ","
                      + std::to_string(findiff) + " result=" + close_sigma.str();
                CHECK_MESSAGE(close_sigma.isclose, msg);
            }
        }
        auto knots = g2f::get_sersic_mix_knots(order);
        double sersicindices[2] = {interpolator.get_sersicindex_min(), interpolator.get_sersicindex_max()};
        g2f::SersicMixValues expected[2] = {knots[0], knots.back()};
        for (size_t i = 0; i < 2; ++i) {
            double sersicindex = sersicindices[i];
            const auto& knots = expected[i];
            CHECK(knots.sersicindex == sersicindex);
            auto interp = interpolator.get_integralsizes(sersicindex);
            for (size_t ord = 0; ord < order; ++ord) {
                CHECK(knots.values[ord].integral == interp[ord].integral);
                CHECK(knots.values[ord].sigma == interp[ord].sigma);
            }
        }
    }
}
