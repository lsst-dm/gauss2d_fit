#include <stdexcept>
#include "componentmixture.h"
#include "psfmodel.h"
#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN

#include "doctest.h"

#include <experimental/iterator>
#include <set>
#include <sstream>

#include "gauss2d/vectorimage.h"

#include "centroidparameters.h"
#include "channel.h"
#include "data.h"
#include "ellipticalcomponent.h"
#include "fractionalintegralmodel.h"
#include "gaussiancomponent.h"
#include "gaussianprior.h"
#include "integralmodel.h"
#include "linearintegralmodel.h"
#include "model.h"
#include "observation.h"
#include "param_defs.h"
#include "param_filter.h"
#include "parameters.h"
#include "sersicmixcomponent.h"
#include "shapeprior.h"
#include "transforms.h"
#include "util.h"

typedef gauss2d::VectorImage<double> Image;
typedef gauss2d::VectorImage<size_t> Indices;
typedef gauss2d::VectorImage<bool> Mask;
typedef g2f::Observation<double, Image, Mask> Observation;
typedef g2f::Data<double, Image, Mask> Data;
typedef g2f::Model<double, Image, Indices, Mask> Model;
typedef std::vector<std::shared_ptr<const g2f::Channel>> Channels;

namespace g2f = gauss2d::fit;

const auto CHANNEL_NONE = g2f::Channel::NONE_PTR();
Channels CHANNELS_NONE = {CHANNEL_NONE};

std::shared_ptr<Data> make_data(Channels channels, const size_t n_x, const size_t n_y,
                                const double sigma_inv_value = 1.) {
    std::vector<std::shared_ptr<const Observation>> observations;
    for (const auto& channel_ptr : channels) {
        auto img = std::make_unique<Image>(n_y, n_x);
        img->fill(0);
        auto err = std::make_unique<Image>(n_y, n_x);
        err->fill(sigma_inv_value);
        auto mask = std::make_unique<Mask>(n_y, n_x);
        mask->fill(1);
        auto observation = std::make_shared<Observation>(std::move(img), std::move(err), std::move(mask),
                                                         *channel_ptr);
        observations.emplace_back(observation);
    }
    return std::make_shared<Data>(observations);
}

g2f::LinearIntegralModel::Data make_integrals(
        const std::vector<std::shared_ptr<const g2f::Channel>>& channels, double value = 1.,
        bool fixed = false) {
    g2f::LinearIntegralModel::Data integrals{};
    for (const auto& channel : channels) {
        auto param = std::make_shared<g2f::IntegralParameter>(
                value, nullptr, g2f::get_transform_default<g2f::Log10Transform>(), nullptr, fixed,
                (*channel).name);
        integrals.emplace_back(*channel, std::move(param));
    }
    return integrals;
}

template <typename t>
std::shared_ptr<t> make_shared_fixed(double value) {
    return std::make_shared<t>(value, nullptr, nullptr, nullptr, true);
}

void verify_model(Model& model, const std::vector<std::shared_ptr<const g2f::Channel>>& channels,
                  g2f::ParamRefs params_free, bool check_outputs_differ = true, bool skip_rho = false,
                  bool print = false, double findiff_frac = 1e-4, double findiff_add = 1e-4,
                  double rtol = 1e-3, double atol = 1e-3) {
    const size_t n_channels = channels.size();
    auto grads = model.compute_loglike_grad(true, false, true);
    g2f::ParamFilter filter_free{false, true, true, true};
    params_free = model.get_parameters_new(&filter_free);
    params_free = g2f::nonconsecutive_unique(params_free);
    CHECK(grads.size() == params_free.size());

    model.setup_evaluators(Model::EvaluatorMode::image);
    auto loglike_img = model.evaluate();
    size_t n_loglike = loglike_img.size();
    for (size_t idx_like = 0; idx_like < n_loglike; ++idx_like) {
        CHECK(loglike_img[idx_like] == 0);
    }

    model.setup_evaluators(Model::EvaluatorMode::loglike);
    auto loglike = model.evaluate();

    for (size_t i = 0; i < 2; ++i) {
        for (unsigned short do_jacobian = false; do_jacobian <= true; do_jacobian++) {
            model.setup_evaluators(
                do_jacobian ? Model::EvaluatorMode::jacobian : Model::EvaluatorMode::loglike_image,
                {}, {}, {}, {}, print
            );

            auto result = model.evaluate();
            for (size_t idx_like; idx_like < n_loglike; ++idx_like) {
                auto close = g2f::isclose(loglike[idx_like], result[idx_like]);
                CHECK_MESSAGE(close.isclose, "loglike[", std::to_string(idx_like), "] !close: ", close.str());
            }

            auto outputs = model.get_outputs();
            if (do_jacobian) {
                auto errors = model.verify_jacobian(findiff_frac, findiff_add, rtol, atol);
                std::string errormsg = "";
                if (errors.size() != 0) {
                    std::stringstream ss;
                    ss << "";
                    std::copy(std::begin(errors), std::end(errors),
                              std::experimental::make_ostream_joiner(ss, "\n"));
                    errormsg = ss.str();
                }
                CHECK(errormsg == "");
            } else {
                CHECK(outputs.size() == n_channels);
                for (size_t idx_channel = 1; idx_channel < n_channels; ++idx_channel) {
                    bool values_equal = outputs[0]->get_value(0, 0) == outputs[idx_channel]->get_value(0, 0);
                    CHECK(values_equal != check_outputs_differ);
                }
            }
        }

        for (auto& param : params_free) {
            double value = param.get().get_value_transformed();
            try {
                param.get().set_value_transformed(value + 0.01);
                param.get().set_value_transformed(value);
            } catch (const std::runtime_error& err) {
                param.get().set_value_transformed(value - 0.01);
                param.get().set_value_transformed(value);
            }
        }
    }

    std::vector<double> values_transformed;
    for (const auto& param : params_free) values_transformed.push_back(param.get().get_value_transformed());
    model.setup_evaluators(Model::EvaluatorMode::image);
    model.evaluate();
    auto outputs = model.get_outputs();
    for (size_t idx = 0; idx < model.size(); ++idx) {
        auto& datum = model.get_data()->at(idx).get().get_image();
        auto output = outputs.at(idx);
        const size_t n_rows = output->get_n_rows();
        const size_t n_cols = output->get_n_cols();
        CHECK(n_rows == datum.get_n_rows());
        CHECK(n_cols == datum.get_n_cols());
        for (size_t row = 0; row < n_rows; ++row) {
            for (size_t col = 0; col < n_cols; ++col) {
                datum.set_value_unchecked(row, col, output->get_value_unchecked(row, col));
            }
        }
    }

    for (unsigned int offset = 0; offset <= 1; ++offset) {
        if (offset) {
            size_t idx_param = 0;
            for (auto& paramref : params_free) {
                auto& param = paramref.get();
                double diff = 1e-4 * (1.0 - 2 * (idx_param % 2));
                diff = g2f::finite_difference_param(param, diff);
                values_transformed[idx_param++] = param.get_value_transformed();
            }
        }
        for (unsigned int include_prior = 0; include_prior <= 1; ++include_prior) {
            for (unsigned int transformed = 0; transformed <= 1; ++transformed) {
                auto hessian = model.compute_hessian(transformed, include_prior);
                size_t idx_param = 0;
                for (const auto& paramref : params_free) {
                    const auto& param = paramref.get();
                    double value_new = param.get_value_transformed();
                    double value_old = values_transformed[idx_param++];
                    CHECK_MESSAGE(g2f::isclose(value_old, value_new, 1e-10, 1e-12).isclose, param.str(),
                                  " value_transformed changed from ", g2f::to_string_float(value_old), " to ",
                                  g2f::to_string_float(value_new),
                                  " with compute_hessian(transformed=", std::to_string(transformed), ")");
                }
                g2f::HessianOptions opt{false, 1e-10, 1e-10};
                auto hessian2 = model.compute_hessian(transformed, include_prior, opt);
                const auto n_cols = hessian->get_n_cols();
                const auto n_rows = hessian->get_n_rows();
                for (size_t row = 0; row < n_rows; ++row) {
                    for (size_t col = 0; col < n_cols; ++col) {
                        auto value = hessian->get_value(row, col);
                        if (row == col) {
                            CHECK_MESSAGE(value <= 0, "hessian[", row, ",", col, "] == 0 (", value, ")");
                        }
                        auto& param_str_row = params_free[row].get();
                        auto& param_str_col = params_free[col].get();
                        std::string name_row = param_str_row.get_name();
                        std::string name_col = param_str_col.get_name();
                        if (!((name_row == g2f::SersicIndexParameter::_name)
                              || (name_row == g2f::ReffXParameter::_name)
                              || (name_row == g2f::ReffYParameter::_name)
                              || (skip_rho && (name_row == g2f::RhoParameter::_name))
                              || (name_col == g2f::SersicIndexParameter::_name)
                              || (name_col == g2f::ReffXParameter::_name)
                              || (name_col == g2f::ReffYParameter::_name)
                              || (skip_rho && (name_col == g2f::RhoParameter::_name)))) {
                            auto value2 = hessian2->get_value(row, col);
                            CHECK_MESSAGE(g2f::isclose(value, value2, 1e-3, 1e-4).isclose, "hessian[", row,
                                          ",", col, "]=(", value, ")!= hessian2 value=(", value2,
                                          ") for params_free[row], params_free[col]=", param_str_row.str(),
                                          param_str_col.str());
                        }
                    }
                    params_free[row].get().set_value_transformed(values_transformed[row]);
                }
            }
        }
    }
}

TEST_CASE("Model") {
    const std::vector<std::shared_ptr<const g2f::Channel>> channels
            = {// TODO: Figure out how this works - auto const conversion?
               g2f::Channel::make("r"), g2f::Channel::make("g"), g2f::Channel::make("b")};

    auto data = make_data(channels, 11, 13);
    g2f::ParamCRefs params_data{};
    data->get_parameters_const(params_data);
    CHECK(params_data.size() == 0);

    Model::PsfModels psfmodels{};
    for (size_t i = 0; i < data->size(); ++i) {
        g2f::Components comps;
        const double integral_factor = 1.1;
        auto integrals = make_integrals(CHANNELS_NONE, integral_factor, true);
        auto model_total = std::make_shared<g2f::LinearIntegralModel>(&integrals);
        g2f::ParamRefs params_integralmodel;
        model_total->get_parameters(params_integralmodel);
        for (const auto& param : params_integralmodel) param.get().set_fixed(true);
        std::shared_ptr<g2f::IntegralModel> last = model_total;

        double integral = 0;
        std::vector<std::pair<double, double>> sizefracs = {{0., 1. / (i + 2.)}, {2., 1.}};
        const size_t idx_sizefrac_last = sizefracs.size() - 1;
        for (size_t idx_comp = 0; idx_comp <= idx_sizefrac_last; ++idx_comp) {
            const auto& sizefrac = sizefracs[idx_comp];
            auto frac = std::make_shared<g2f::ProperFractionParameter>(sizefrac.second);
            CHECK(frac->get_value() == sizefrac.second);
            // Would set this condition if we wanted the other fractions free
            // if(sizefrac.second == 1)
            frac->set_fixed(true);
            g2f::FractionalIntegralModel::Data data{{g2f::Channel::NONE(), frac}};
            auto model_frac = g2f::FractionalIntegralModel::make(data, last, idx_comp == idx_sizefrac_last);

            double integral_comp = model_frac->get_integral(g2f::Channel::NONE());
            if (idx_comp == 0) {
                CHECK(integral_comp == integral_factor * sizefrac.second);
            }

            integral += integral_comp;

            auto comp = std::make_unique<g2f::GaussianComponent>(
                    std::make_shared<g2f::GaussianParametricEllipse>(sizefrac.first, sizefrac.first, 0),
                    nullptr, model_frac);
            auto gaussians = comp->get_gaussians(g2f::Channel::NONE());
            CHECK(gaussians->size() == 1);
            CHECK(gaussians->at(0).get_integral_value() == integral_comp);

            g2f::ParamRefs params_comp;
            comp->get_parameters(params_comp);
            for (const auto& param : params_comp) param.get().set_fixed(true);
            comps.emplace_back(std::move(comp));
            last = model_frac;
        }
        CHECK(std::abs(integral - integral_factor) < 1e-12);

        auto psfmodel = std::make_shared<g2f::PsfModel>(comps);

        auto params_psf = psfmodel->get_parameters_const_new();
        // 2 comps x (6 gauss + 1 frac) (last constant unity frac omitted)
        CHECK(params_psf.size() == 14);
        for (const auto& param : params_psf) {
            CHECK(param.get().get_fixed() == true);
        }

        const auto gaussians = psfmodel->get_gaussians();
        CHECK(gaussians->size() == 2);

        size_t g = 0;
        for (const auto& gauss : *gaussians) {
            std::string msg = "obs[" + std::to_string(i) + "," + std::to_string(g) + "]";
            CHECK_MESSAGE(gauss->get_integral_value() > 0, msg);
            g++;
        }
        // std::cerr << psfmodel->str() << std::endl;
        // std::cerr << gaussians->str() << std::endl;

        psfmodels.push_back(std::move(psfmodel));
    }

    const bool free_sersicindex = true;

    auto limits_axrat_logit = std::make_shared<parameters::Limits<double>>(-1e-10, 1 + 1e-10);
    auto transform_axrat = std::make_shared<g2f::LogitLimitedTransform>(limits_axrat_logit);

    Model::Sources sources{};
    std::vector<std::shared_ptr<g2f::Prior>> priors = {};
    for (size_t i = 0; i < 2; ++i) {
        std::vector<std::shared_ptr<g2f::Component>> comps;
        for (size_t c = 0; c < 2; ++c) {
            auto centroids = std::make_shared<g2f::CentroidParameters>(c + i + 1, c + i + 1.5);
            priors.emplace_back(std::make_shared<g2f::GaussianPrior>(centroids->get_x_param_ptr(),
                                                                     centroids->get_x(), 1.0, false));
            priors.emplace_back(std::make_shared<g2f::GaussianPrior>(centroids->get_y_param_ptr(),
                                                                     centroids->get_y(), 0.5, false));
            auto integrals = make_integrals(channels, c + 1);
            auto integralmodel = std::make_shared<g2f::LinearIntegralModel>(&integrals);
            std::shared_ptr<g2f::Component> comp;
            std::shared_ptr<g2f::ParametricEllipse> ellipse;
            if (i == 0) {
                auto ellipse_g = std::make_shared<g2f::GaussianParametricEllipse>(c + 0.5, c + 1.5, 0);
                ellipse = ellipse_g;
                comp = std::make_shared<g2f::GaussianComponent>(ellipse_g, centroids, integralmodel);
            } else {
                auto sersic_n = std::make_shared<g2f::SersicMixComponentIndexParameter>(0.5 + 3.5 * c);
                sersic_n->set_free(free_sersicindex);
                auto reff_x = std::make_shared<g2f::ReffXParameter>(
                        c + 0.5, nullptr, g2f::get_transform_default<g2f::Log10Transform>());
                auto reff_y = std::make_shared<g2f::ReffYParameter>(
                        c + 1.5, nullptr, g2f::get_transform_default<g2f::Log10Transform>());
                auto ellipse_s = std::make_shared<g2f::SersicParametricEllipse>(reff_x, reff_y);
                ellipse = ellipse_s;
                comp = std::make_shared<g2f::SersicMixComponent>(ellipse_s, centroids, integralmodel,
                                                                 sersic_n);
            }
            auto ell_g2 = gauss2d::Ellipse(ellipse->get_size_x(), ellipse->get_size_y(), ellipse->get_rho());
            auto ell_maj = gauss2d::EllipseMajor(ell_g2);
            double axrat = ell_maj.get_axrat();
            double size_ell = sqrt(std::pow(ell_maj.get_r_major(), 2) * axrat
                                   + std::pow(g2f::ShapePriorOptions::size_maj_floor_default, 2));
            axrat = sqrt(axrat * axrat + std::pow(g2f::ShapePriorOptions::axrat_floor_default, 2));

            auto prior_size = std::make_shared<g2f::ParametricGaussian1D>(
                    std::make_shared<g2f::MeanParameter>(size_ell, nullptr,
                                                         g2f::get_transform_default<g2f::Log10Transform>()),
                    std::make_shared<g2f::StdDevParameter>(0.5));
            auto prior_axrat = std::make_shared<g2f::ParametricGaussian1D>(
                    std::make_shared<g2f::MeanParameter>(axrat, nullptr, transform_axrat),
                    std::make_shared<g2f::StdDevParameter>(1.0));
            priors.emplace_back(std::make_shared<g2f::ShapePrior>(ellipse, prior_size, prior_axrat));
            comps.emplace_back(comp);
        }
        auto source = std::make_shared<g2f::Source>(comps);
        sources.push_back(source);
    }
    auto params_src = sources[0]->get_parameters_const_new();
    auto model = std::make_shared<Model>(data, psfmodels, sources, priors);
    CHECK(model->str() != "");
    CHECK(model->get_priors().size() == priors.size());

    auto params = model->get_parameters_new();
    // 2 sources x (2 comps x (3 integral, 2 centroid, 3 ellipse)) = 32
    // + 1 source x (2 comps x 1 sersic_n) = 34
    const size_t n_params_src = 34;
    // PSF: 3 observations x (2 comp x (1 integral, n frac, 2 centroid, 3 ellipse)) = 45
    // (each comp after second has an extra frac per channel)
    const size_t n_params_psf = 42;
    CHECK(params.size() == (n_params_src + n_params_psf));

    // PSF: 3 observations x (1 integral + 1 frac + 2 comp x (2 centroid, 3 ellipse)) = 36
    std::set<g2f::ParamBaseCRef> paramset(params.cbegin(), params.cend());
    const size_t n_params_psf_uniq = 36;
    CHECK(paramset.size() == n_params_src + n_params_psf_uniq);

    g2f::ParamFilter filter_free{false, true, true, true};
    params = model->get_parameters_new(&filter_free);
    params = g2f::nonconsecutive_unique(params);
    CHECK(params.size() == n_params_src);

    const auto& channel = *channels[0];

    CHECK(model->get_n_gaussians(channel) == 10);
    auto gaussians = model->get_gaussians(channel);
    const size_t n_gauss_src = gaussians->size();
    // 2 comps x 1 (gauss) + 2 comps x 4 (sersic)
    CHECK(n_gauss_src == 10);

    const auto psfcomps = model->get_psfmodels().at(0)->get_gaussians();
    const size_t n_gauss_psf = psfcomps->size();
    CHECK(n_gauss_psf == 2);
    for (size_t p = 0; p < n_gauss_psf; ++p) CHECK(psfcomps->at(p).get_integral_value() > 0);

    auto map_extra = std::make_shared<g2f::ExtraParamMap>();
    auto map_grad = std::make_shared<g2f::GradParamMap>();

    auto factors_extra = std::make_shared<g2f::ExtraParamFactors>();
    auto factors_grad = std::make_shared<g2f::GradParamFactors>();

    const size_t n_gauss = n_gauss_src * n_gauss_psf;

    map_extra->reserve(n_gauss);
    map_grad->reserve(n_gauss);
    factors_extra->reserve(n_gauss);
    factors_grad->reserve(n_gauss);

    // Check that filtering by channel works
    g2f::ParamCRefs params_src_free_const;
    g2f::ParamRefs params_src_free;
    g2f::ParamFilter filter{false, true, true, true, channel};
    g2f::ParameterMap offsets{};

    for (const auto& source : model->get_sources()) {
        source->get_parameters_const(params_src_free_const, &filter);
        source->get_parameters(params_src_free, &filter);
        for (size_t i_psf = 0; i_psf < n_gauss_psf; ++i_psf) {
            source->add_grad_param_map(channel, *map_grad, offsets);
            source->add_extra_param_map(channel, *map_extra, *map_grad, offsets);

            source->add_grad_param_factors(channel, *factors_grad);
            source->add_extra_param_factors(channel, *factors_extra);
        }
    }
    auto n_free = params_src_free.size();
    // 2 sources x (2 comps x (1 integral, 2 centroid, 3 ellipse)) = 24
    // + 1 source x (2 comps x 1 sersicindex) = 26
    CHECK(n_free == 24 + 2 * free_sersicindex);

    CHECK(offsets.size() == n_free);

    CHECK(map_extra->size() == n_gauss);
    CHECK(map_grad->size() == n_gauss);
    CHECK(factors_extra->size() == n_gauss);
    CHECK(factors_grad->size() == n_gauss);

    // TODO: This should pass without setting skip_rho
    verify_model(*model, channels, params_src_free, true, true);

    // Add a fractional source, which should work even if it's not useful
    auto centroids = std::make_shared<g2f::CentroidParameters>(5.5, 6.5);
    g2f::FractionalIntegralModel::Data data_frac1 = {};
    g2f::FractionalIntegralModel::Data data_frac2 = {};
    double frac_value = 0.3;
    for (const auto& channel : channels) {
        auto frac1 = std::make_shared<g2f::ProperFractionParameter>(frac_value);
        auto frac2 = std::make_shared<g2f::ProperFractionParameter>(1.0, nullptr, nullptr, nullptr, true);
        data_frac1.emplace_back(*channel, frac1);
        data_frac2.emplace_back(*channel, frac2);
        frac_value += 0.15;
    }

    auto integrals_frac = make_integrals(channels, 0.3);
    auto integralmodel_frac = std::make_shared<g2f::LinearIntegralModel>(&integrals_frac);

    auto fracmodel1 = g2f::FractionalIntegralModel::make(data_frac1, integralmodel_frac, false);
    g2f::Components comps_src3 = {
            std::make_shared<g2f::GaussianComponent>(
                    std::make_shared<g2f::GaussianParametricEllipse>(2.3, 4.4, 0.15), centroids, fracmodel1),
            std::make_shared<g2f::GaussianComponent>(
                    std::make_shared<g2f::GaussianParametricEllipse>(2.6, 3.8, -0.33), centroids,
                    g2f::FractionalIntegralModel::make(data_frac2, fracmodel1, true))};
    auto source = std::make_shared<g2f::Source>(comps_src3);
    sources.clear();
    sources.push_back(source);

    auto psfmodels_ref = model->get_psfmodels();
    for (size_t idx = 0; idx < psfmodels_ref.size(); ++idx) {
        psfmodels[idx] = psfmodels_ref[idx];
    }
    priors.clear();
    auto model2 = std::make_shared<Model>(data, psfmodels, sources, priors);
    auto params_src_free2 = model2->get_parameters_new(&filter_free);
    params_src_free2 = g2f::nonconsecutive_unique<g2f::ParamBaseRef>(params_src_free2);

    // TODO: enable when DM-40674 is fixed
    // verify_model(*model2, channels, params_src_free2, true, true);
}

TEST_CASE("Model PSF") {
    const Channels channels = CHANNELS_NONE;
    auto data = make_data(channels, 11, 13);
    std::vector<std::shared_ptr<g2f::Component>> comps = {};
    auto integrals = make_integrals(CHANNELS_NONE, 7.0, true);
    auto model_total = std::make_shared<g2f::LinearIntegralModel>(&integrals);
    std::shared_ptr<g2f::IntegralModel> last = model_total;

    auto cens = std::make_shared<g2f::CentroidParameters>(std::make_shared<g2f::CentroidXParameter>(0),
                                                          std::make_shared<g2f::CentroidYParameter>(0));

    for (const auto& sizefrac : {std::pair{1.5, 1.0 - 1e-15}, {2.5, 1.0}}) {
        const auto is_last = sizefrac.second == 1;
        auto frac = std::make_shared<g2f::ProperFractionParameter>(sizefrac.second);
        // Would set this condition if we wanted the other fractions free
        if (is_last) frac->set_fixed(true);
        g2f::FractionalIntegralModel::Data data{{g2f::Channel::NONE(), frac}};
        const double& size = sizefrac.first;

        last = g2f::FractionalIntegralModel::make(data, last, is_last);

        comps.emplace_back(std::make_shared<g2f::GaussianComponent>(
                std::make_shared<g2f::GaussianParametricEllipse>(std::make_shared<g2f::SigmaXParameter>(size),
                                                                 std::make_shared<g2f::SigmaYParameter>(size),
                                                                 std::make_shared<g2f::RhoParameter>(0)),
                cens, last));
    }
    g2f::Components psfcomps = {g2f::GaussianComponent::make_uniq_default_gaussians({0.})};
    Model::PsfModels psfmodels({std::make_shared<g2f::PsfModel>(psfcomps)});
    Model::Sources sources({std::make_shared<g2f::Source>(comps)});

    std::vector<std::shared_ptr<g2f::Prior>> priors = {};
    auto model = std::make_shared<Model>(data, psfmodels, sources, priors);
    g2f::ParamRefs params_free;
    g2f::ParamFilter filter{false, true, true, true, g2f::Channel::NONE()};
    g2f::ParameterMap offsets{};

    model->get_parameters(params_free, &filter);
    params_free = g2f::nonconsecutive_unique(params_free);
    // 2 cens + 1 fluxfrac + 2*(2sigmas + rho) = 9
    CHECK(params_free.size() == 9);
    verify_model(*model, channels, params_free, true, false, true);
}

// This has some overlap with the Model test case but it's not really problematic
TEST_CASE("Model with priors") {
    const size_t n_x = 33, n_y = 33;
    auto data = make_data(CHANNELS_NONE, n_x, n_y);
    auto centroid_psf = std::make_shared<g2f::CentroidParameters>(
            make_shared_fixed<g2f::CentroidXParameter>(0), make_shared_fixed<g2f::CentroidYParameter>(0));
    g2f::LinearIntegralModel::Data model_psf_ref_data
            = {{*CHANNEL_NONE, make_shared_fixed<g2f::IntegralParameter>(1.0)}};
    auto model_psf_ref = std::make_shared<g2f::LinearIntegralModel>(&model_psf_ref_data);

    g2f::FractionalIntegralModel::Data model_psf_flux_first_data
            = {{*CHANNEL_NONE, make_shared_fixed<g2f::ProperFractionParameter>(0.806268)}};
    auto model_psf_flux_first
            = g2f::FractionalIntegralModel::make(model_psf_flux_first_data, model_psf_ref, false);

    g2f::FractionalIntegralModel::Data model_psf_flux_second_data
            = {{*CHANNEL_NONE, make_shared_fixed<g2f::ProperFractionParameter>(1.0)}};
    auto model_psf_flux_second
            = g2f::FractionalIntegralModel::make(model_psf_flux_second_data, model_psf_flux_first, true);

    g2f::Components comps_psf = {std::make_shared<g2f::GaussianComponent>(
                                         std::make_shared<g2f::GaussianParametricEllipse>(
                                                 make_shared_fixed<g2f::SigmaXParameter>(1.648791),
                                                 make_shared_fixed<g2f::SigmaYParameter>(1.619646),
                                                 make_shared_fixed<g2f::RhoParameter>(-0.009130)),
                                         centroid_psf, model_psf_flux_first),
                                 std::make_shared<g2f::GaussianComponent>(
                                         std::make_shared<g2f::GaussianParametricEllipse>(
                                                 make_shared_fixed<g2f::SigmaXParameter>(3.325370),
                                                 make_shared_fixed<g2f::SigmaYParameter>(3.346987),
                                                 make_shared_fixed<g2f::RhoParameter>(-0.008666)),
                                         centroid_psf, model_psf_flux_first)};
    auto psfmodel = std::make_shared<g2f::PsfModel>(comps_psf);
    Model::PsfModels psfmodels = {psfmodel};

    g2f::LinearIntegralModel::Data model_src_flux_data
            = {{*CHANNEL_NONE, std::make_shared<g2f::IntegralParameter>(0.367757)}};

    auto ellipse_src = std::make_shared<g2f::SersicParametricEllipse>(
            std::make_shared<g2f::ReffXParameter>(1.475109), std::make_shared<g2f::ReffYParameter>(1.469544),
            std::make_shared<g2f::RhoParameter>(-0.020842));

    auto ell_g2
            = gauss2d::Ellipse(ellipse_src->get_size_x(), ellipse_src->get_size_y(), ellipse_src->get_rho());
    auto ell_maj = gauss2d::EllipseMajor(ell_g2);
    double axrat = ell_maj.get_axrat();
    double size_ell = sqrt(std::pow(ell_maj.get_r_major(), 2) * axrat);
    axrat = sqrt(axrat * axrat + std::pow(g2f::ShapePriorOptions::axrat_floor_default, 2));

    g2f::Components comps_src = {std::make_shared<g2f::SersicMixComponent>(
            ellipse_src,
            std::make_shared<g2f::CentroidParameters>(make_shared_fixed<g2f::CentroidXParameter>(16.586227),
                                                      make_shared_fixed<g2f::CentroidYParameter>(16.695350)),
            std::make_shared<g2f::LinearIntegralModel>(&model_src_flux_data),
            make_shared_fixed<g2f::SersicMixComponentIndexParameter>(1.0))};

    Model::Sources sources = {std::make_shared<g2f::Source>(comps_src)};

    Model::Priors priors = {std::make_shared<g2f::ShapePrior>(
            ellipse_src,
            std::make_shared<g2f::ParametricGaussian1D>(std::make_shared<g2f::MeanParameter>(size_ell),
                                                        std::make_shared<g2f::StdDevParameter>(0.3)),
            std::make_shared<g2f::ParametricGaussian1D>(std::make_shared<g2f::MeanParameter>(axrat),
                                                        std::make_shared<g2f::StdDevParameter>(0.2)))};

    auto model = std::make_shared<Model>(data, psfmodels, sources, priors);

    g2f::ParamFilter filter{false, true, true, true, g2f::Channel::NONE()};
    auto params_free = model->get_parameters_new(&filter);
    params_free = g2f::nonconsecutive_unique(params_free);

    verify_model(*model, CHANNELS_NONE, params_free, true, false, true, 1e-5, 1e-5);

    // Test repeat evaluation with and without forcing and specifiying outputs
    model->setup_evaluators(Model::EvaluatorMode::jacobian);
    model->evaluate();

    model->setup_evaluators(Model::EvaluatorMode::loglike_image);
    model->evaluate();
    model->setup_evaluators(Model::EvaluatorMode::jacobian, {}, {}, {}, nullptr, false);
    model->evaluate();
    std::vector<std::vector<std::shared_ptr<Image>>> outputs = {{std::make_shared<Image>(n_y, n_x)}};
    model->setup_evaluators(Model::EvaluatorMode::loglike_image, outputs, {}, {}, nullptr, false);
    model->evaluate();
    // outputs is too short
    CHECK_THROWS_AS(model->setup_evaluators(Model::EvaluatorMode::jacobian, outputs, {}, {}, nullptr, true),
                    std::invalid_argument);
    // test that failed setup causes evaluate to throw until set up successfully
    CHECK_THROWS_AS(model->evaluate(), std::runtime_error);
    outputs[0].resize(params_free.size() + 1, outputs[0][0]);
    model->setup_evaluators(Model::EvaluatorMode::jacobian, outputs, {}, {}, nullptr, true);
    model->evaluate();
}
