#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN

#include "doctest.h"

#include <cmath>
#include <memory>

#include "lsst/gauss2d/fit/gaussianparametricellipse.h"
#include "lsst/gauss2d/fit/parameters.h"
#include "lsst/gauss2d/fit/parametricgaussian1d.h"
#include "lsst/gauss2d/fit/shapeprior.h"
#include "lsst/gauss2d/fit/transforms.h"
#include "lsst/gauss2d/fit/util.h"

namespace g2d = lsst::gauss2d;
namespace g2f = lsst::gauss2d::fit;

TEST_CASE("ShapePrior") {
    CHECK_THROWS_AS(g2f::ShapePrior(nullptr), std::invalid_argument);
    auto options = std::make_shared<g2f::ShapePriorOptions>();
    CHECK_GT(options->repr().size(), 0);
    CHECK_GT(options->repr(true).size(), 0);
    CHECK_GT(options->str().size(), 0);
    double mean_size = 16.6;
    double stddev_size = 3.4;
    double mean_axrat = 0.5;
    double stddev_axrat = 0.1;

    auto ellipse = std::make_shared<g2f::GaussianParametricEllipse>();
    auto ellipse_mean = std::make_shared<g2d::Ellipse>(g2d::EllipseMajor(mean_size, mean_axrat, 0));
    ellipse->set_xyr(ellipse_mean->get_xyr());

    CHECK_GT(g2f::ShapePrior(ellipse).str().size(), 0);
    CHECK_GT(g2f::ShapePrior(ellipse, nullptr, nullptr, options).repr().size(), 0);

    auto limits_axrat_logit = std::make_shared<parameters::Limits<double>>(-1e-10, 1 + 1e-10);
    auto transform_axrat = std::make_shared<g2f::LogitLimitedTransform>(limits_axrat_logit);

    auto prior_size_params = std::make_shared<g2f::ParametricGaussian1D>(
            std::make_shared<g2f::MeanParameterD>(mean_size, nullptr,
                                                  g2f::get_transform_default<g2f::Log10Transform>()),
            std::make_shared<g2f::StdDevParameterD>(stddev_size));
    auto prior_axrat_params = std::make_shared<g2f::ParametricGaussian1D>(
            std::make_shared<g2f::MeanParameterD>(mean_axrat, nullptr, transform_axrat),
            std::make_shared<g2f::StdDevParameterD>(stddev_axrat));

    auto prior_size = g2f::ShapePrior(ellipse, prior_size_params, nullptr, options);
    double loglike_size = prior_size.evaluate().loglike;
    CHECK_NE(loglike_size, 0);

    double delta_param = 1e-6;
    double mean_size_trans = prior_size_params->get_mean_parameter().get_value_transformed();
    prior_size_params->get_mean_parameter().set_value_transformed(mean_size_trans + delta_param / 2.);
    double loglike_plus = prior_size.evaluate().loglike;
    CHECK_LT(loglike_plus, loglike_size);
    prior_size_params->get_mean_parameter().set_value_transformed(mean_size_trans - delta_param / 2.);
    double loglike_minus = prior_size.evaluate().loglike;
    // TODO: Rethink this for product size, if/when implemented
    auto isclose_size = g2f::isclose(loglike_plus, loglike_minus, 1e-5, 1e-7);
    CHECK_MESSAGE(isclose_size.isclose, "is_close=", isclose_size.str());
    prior_size_params->get_mean_parameter().set_value(mean_size);

    auto prior_axrat = g2f::ShapePrior(ellipse, nullptr, prior_axrat_params, options);
    double loglike_axrat = prior_axrat.evaluate().loglike;
    CHECK_NE(loglike_axrat, 0);

    delta_param = 1e-4;
    double mean_axrat_trans = prior_axrat_params->get_mean_parameter().get_value_transformed();
    prior_axrat_params->get_mean_parameter().set_value_transformed(mean_axrat_trans + delta_param / 2.);
    loglike_plus = prior_axrat.evaluate().loglike;
    CHECK_LT(loglike_plus, loglike_axrat);
    prior_axrat_params->get_mean_parameter().set_value_transformed(mean_axrat_trans - delta_param / 2.);
    loglike_minus = prior_axrat.evaluate().loglike;
    auto isclose_axrat = g2f::isclose(loglike_plus, loglike_minus, 1e-5, 1e-7);
    CHECK_MESSAGE(isclose_axrat.isclose, "is_close=", isclose_axrat.str());
    prior_axrat_params->get_mean_parameter().set_value_transformed(mean_axrat_trans);

    auto prior = g2f::ShapePrior(ellipse, prior_size_params, prior_axrat_params, options);

    CHECK_EQ(g2f::isclose(prior.evaluate().loglike, loglike_axrat + loglike_size, 1e-8, 1e-10).isclose, true);
    auto terms = prior.get_loglike_const_terms();
    CHECK_EQ(terms.size(), 2);
    CHECK_GT(prior.repr().size(), 0);
    CHECK_GT(prior.repr(true).size(), 0);
    CHECK_GT(prior.str().size(), 0);

    std::map<g2f::ParamBaseRef, std::vector<double>> param_values_test
            = {{ellipse->get_sigma_x_param(), {0.5, 2.4}},
               {ellipse->get_sigma_y_param(), {0.13, 1.0}},
               {ellipse->get_rho_param(), {-0.83, 0., 0.18}}};

    delta_param = 1e-5;

    for (auto& [paramref, values] : param_values_test) {
        for (double x : values) {
            auto& param = paramref.get();
            param.set_value_transformed(x - delta_param / 2.);
            auto eval_minus = prior.evaluate();
            param.set_value_transformed(x + delta_param / 2.);
            auto eval_plus = prior.evaluate();
            const double dll_dx_findiff = (eval_plus.loglike - eval_minus.loglike) / delta_param;
            param.set_value_transformed(x);
            auto eval = prior.evaluate(true);
            const double dll_dx_est = eval.compute_dloglike_dx(param);

            auto result = g2f::isclose(dll_dx_findiff, dll_dx_est, 1e-6, 1e-5);
            CHECK_MESSAGE(result.isclose, param.str(), "is_close=", result.str());
            for (size_t idx_resid = 0; idx_resid < eval.residuals.size(); idx_resid++) {
                double jacob = eval.jacobians[paramref][idx_resid];
                const double dresid_dx_findiff
                        = (eval_plus.residuals[idx_resid] - eval_minus.residuals[idx_resid]) / delta_param;
                auto result_resid = g2f::isclose(dresid_dx_findiff, jacob, 1e-6, 1e-5);
                CHECK_MESSAGE(result_resid.isclose, param.str(), "is_close=", result_resid.str());
            }
        }
    }
}
