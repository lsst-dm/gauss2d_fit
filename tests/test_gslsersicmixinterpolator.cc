#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN

#include "doctest.h"

#ifdef LSST_GAUSS2D_FIT_HAS_GSL

#include "lsst/gauss2d/fit/gslsersicmixinterpolator.h"
#include "lsst/gauss2d/fit/util.h"

namespace g2f = lsst::gauss2d::fit;

TEST_CASE("GSLSersicMixInterpolator") {
    auto interpolator_default = g2f::GSLSersicMixInterpolator();
    CHECK_EQ(interpolator_default.get_order(), g2f::SERSICMIX_ORDER_DEFAULT);
    CHECK_EQ(interpolator_default.get_interptype(), g2f::GSLInterpolator::INTERPTYPE_DEFAULT);
    CHECK_GT(interpolator_default.str().size(), 0);

    std::vector<unsigned short> orders = {4, 8};
    for (const auto order : orders) {
        auto interpolator = g2f::GSLSersicMixInterpolator(order);
        CHECK_GT(interpolator.get_sersicindex_max(), interpolator.get_sersicindex_min());
        CHECK_THROWS_AS(interpolator.get_integralsizes(0.499), std::invalid_argument);
        CHECK_THROWS_AS(interpolator.get_integralsizes_derivs(0.499), std::invalid_argument);
        double too_big = interpolator.get_sersicindex_max() + 1e-10;
        CHECK_THROWS_AS(interpolator.get_integralsizes(too_big), std::invalid_argument);
        CHECK_THROWS_AS(interpolator.get_integralsizes_derivs(too_big), std::invalid_argument);

        double diff = interpolator.get_sersicindex_max() - interpolator.get_sersicindex_min();

        for (const double factor : {0., 0.5, 1.}) {
            double delta = 1e-8;
            const double sersicindex = interpolator.get_sersicindex_min() + diff * factor;
            auto result = interpolator.get_integralsizes(sersicindex);
            CHECK_EQ(result.size(), order);
            for (size_t ord = 0; ord < order; ++ord) {
                CHECK_GE(result[ord].integral, 0);
                CHECK_GE(result[ord].sigma, 0);
            }

            double factor_new = factor + delta;
            if (!(factor_new <= 1)) {
                delta = -delta;
                factor_new = factor + delta;
            }
            const double sersicindex_new = interpolator.get_sersicindex_min() + diff * factor_new;
            const double dsersicindex = sersicindex_new - sersicindex;
            auto result_new = interpolator.get_integralsizes(sersicindex_new);
            // This is finite differencing at (x + delta/2) +/- delta/2, which is more accurate
            // (the non-GSL linear interpolator should do this too)
            auto result_derivs = interpolator.get_integralsizes_derivs((sersicindex + sersicindex_new) / 2.);
            for (size_t ord = 0; ord < order; ++ord) {
                CHECK_GE(result_new[ord].integral, 0);
                CHECK_GE(result_new[ord].sigma, 0);
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
        auto knots = lsst::gauss2d::fit::get_sersic_mix_knots(order);
        double sersicindices[2] = {interpolator.get_sersicindex_min(), interpolator.get_sersicindex_max()};
        g2f::SersicMixValues expected[2] = {knots[0], knots.back()};
        for (size_t i = 0; i < 2; ++i) {
            double sersicindex = sersicindices[i];
            const auto& knots = expected[i];
            CHECK_EQ(knots.sersicindex, sersicindex);
            auto interp = interpolator.get_integralsizes(sersicindex);
            for (size_t ord = 0; ord < order; ++ord) {
                if (interpolator.correct_final_integral && ((ord + 1) == order)) {
                    double correction = interpolator.get_final_correction();
                    CHECK_EQ(g2f::isclose(knots.values[ord].integral, correction * interp[ord].integral)
                                     .isclose,
                             true);
                } else {
                    CHECK_EQ(knots.values[ord].integral, interp[ord].integral);
                }
                CHECK_EQ(knots.values[ord].sigma, interp[ord].sigma);
            }
        }
    }
}

#endif