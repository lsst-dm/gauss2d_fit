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
#include "integralmodel.h"
#include "linearintegralmodel.h"
#include "model.h"
#include "observation.h"
#include "param_defs.h"
#include "param_filter.h"
#include "parameters.h"
#include "sersicmixcomponent.h"
#include "util.h"

typedef gauss2d::VectorImage<double> Image;
typedef gauss2d::VectorImage<size_t> Indices;
typedef gauss2d::VectorImage<bool> Mask;
typedef g2f::Observation<double, Image, Mask> Observation;
typedef g2f::Data<double, Image, Mask> Data;
typedef g2f::Model<double, Image, Indices, Mask> Model;
typedef std::vector<std::shared_ptr<const g2f::Channel>> Channels;

namespace g2f = gauss2d::fit;

Channels CHANNELS_NONE = {g2f::Channel::NONE_PTR()};

std::shared_ptr<Data> make_data(
    Channels channels,
    const size_t n_x,
    const size_t n_y,
    const double sigma_inv_value=1.
) {
    std::vector<std::shared_ptr<const Observation>> observations;
    for(const auto & channel_ptr : channels)
    {
        auto err = std::make_unique<Image>(n_y, n_x);
        *err  += sigma_inv_value;
        auto observation = std::make_shared<Observation>(
            std::make_unique<Image>(n_y, n_x),
            std::move(err),
            std::make_unique<Mask>(n_y, n_x),
            *channel_ptr
        );
        observations.emplace_back(observation);
    }
    return std::make_shared<Data>(observations);
}

g2f::LinearIntegralModel::Data make_integrals(
    const std::vector<std::shared_ptr<const g2f::Channel>> & channels,
    double value = 1.,
    bool fixed = false
) {
    g2f::LinearIntegralModel::Data integrals{};
    for(const auto & channel: channels) {
        auto param = std::make_shared<g2f::IntegralParameter>(value);
        if(fixed) param->set_fixed(true);
        integrals[*channel] = std::move(param);
    }
    return integrals;
}

void verify_model(
    Model & model,
    const std::vector<std::shared_ptr<const g2f::Channel>> & channels,
    g2f::ParamRefs params_free,
    bool check_outputs_differ=true,
    bool print=false
) {
    const size_t n_channels = channels.size();
    for(size_t i = 0; i < 2; ++i)
    {
        for(unsigned short do_jacobian=false; do_jacobian <= true; do_jacobian++)
        {            
            model.setup_evaluators(
                do_jacobian ? Model::EvaluatorMode::jacobian : Model::EvaluatorMode::image, {}, {}, print
            );

            std::vector<double> result = model.evaluate();

            auto outputs = model.get_outputs();
            if(do_jacobian) {
                auto errors = model.verify_jacobian();
                std::string errormsg = "";
                if(errors.size() != 0) {
                    std::stringstream ss;
                    ss << "";
                    std::copy(
                        std::begin(errors),
                        std::end(errors),
                        std::experimental::make_ostream_joiner(ss, "\n")
                    );
                    errormsg = ss.str();
                }
                CHECK(errormsg == "");
            } else {
                CHECK(outputs.size() == n_channels);
                CHECK(result[0] == 0);
                for(size_t idx_channel = 1; idx_channel < n_channels; ++idx_channel) {
                    CHECK(result[idx_channel] == 0);
                    bool values_equal = outputs[0]->get_value(0, 0) == outputs[1]->get_value(0, 0);
                    CHECK(values_equal != check_outputs_differ);
                }
            }
        }

        for(auto & param : params_free) {
            double value = param.get().get_value_transformed();
            try {
                param.get().set_value_transformed(value + 0.01);
            } catch(const std::runtime_error & err) {
                param.get().set_value_transformed(value - 0.01);
            }
        }
    }

}

TEST_CASE("Model") {
    const std::vector<std::shared_ptr<const g2f::Channel>> channels = {
        // TODO: Figure out how this works - auto const conversion?
        g2f::Channel::make("r"),
        g2f::Channel::make("g"),
        g2f::Channel::make("b")
    };

    auto data = make_data(channels, 11, 13);
    g2f::ParamCRefs params_data{};
    data->get_parameters_const(params_data);
    CHECK(params_data.size() == 0);

    Model::PsfModels psfmodels{};
    for(size_t i = 0; i < data->size(); ++i) {
        g2f::Components comps;
        auto integrals = make_integrals(CHANNELS_NONE, 1.0, true);
        auto model_total = std::make_shared<g2f::LinearIntegralModel>(&integrals);
        g2f::ParamRefs params_integralmodel;
        model_total->get_parameters(params_integralmodel);
        for(const auto & param : params_integralmodel) param.get().set_fixed(true);
        std::shared_ptr<g2f::IntegralModel> last = model_total;

        double integral = 0;
        std::vector<std::pair<double, double>> sizefracs = {{0., 1./(i + 2.)}, {2., 1.}};
        const size_t idx_sizefrac_last = sizefracs.size() - 1;
        for(size_t idx_comp = 0; idx_comp <= idx_sizefrac_last; ++idx_comp)
        {
            const auto & sizefrac = sizefracs[idx_comp];
            auto frac = std::make_shared<g2f::ProperFractionParameter>(sizefrac.second);
            CHECK(frac->get_value() == sizefrac.second);
            // Would set this condition if we wanted the other fractions free
            // if(sizefrac.second == 1)
            frac->set_fixed(true);
            g2f::FractionalIntegralModel::Data data {{g2f::Channel::NONE(), frac}};
            auto model_frac = g2f::FractionalIntegralModel::make(data, last, idx_comp == idx_sizefrac_last);

            double integral_comp = model_frac->get_integral(g2f::Channel::NONE());
            if(idx_comp == 0) {
                CHECK(integral_comp == sizefrac.second);
            }
            
            integral += integral_comp;
            
            auto comp = std::make_unique<g2f::GaussianComponent>(
                std::make_shared<g2f::GaussianParametricEllipse>(sizefrac.first, sizefrac.first, 0),
                nullptr,
                model_frac
            );
            auto gaussians = comp->get_gaussians(g2f::Channel::NONE());
            CHECK(gaussians->size() == 1);
            CHECK(gaussians->at(0).get_integral_value() == integral_comp);

            g2f::ParamRefs params_comp;
            comp->get_parameters(params_comp);
            for(const auto & param : params_comp) param.get().set_fixed(true);
            comps.emplace_back(std::move(comp));
            last = model_frac;
        }
        CHECK(std::abs(integral - 1) < 1e-12);

        auto psfmodel = std::make_shared<g2f::PsfModel>(comps);

        g2f::ParamCRefs params_psf{};
        psfmodel->get_parameters_const(params_psf);
        // 2 comps x (6 gauss + 1 frac) (last constant unity frac omitted)
        CHECK(params_psf.size() == 14);
        for(const auto & param : params_psf) CHECK(param.get().get_fixed() == true);

        const auto gaussians = psfmodel->get_gaussians();
        CHECK(gaussians->size() == 2);

        size_t g = 0;
        for(const auto & gauss : *gaussians) {
            std::string msg = "obs[" + std::to_string(i) + "," + std::to_string(g) + "]";
            CHECK_MESSAGE(gauss->get_integral_value() > 0, msg);
            g++;
        }
        //std::cerr << psfmodel->str() << std::endl;
        //std::cerr << gaussians->str() << std::endl;

        psfmodels.push_back(std::move(psfmodel));
    }

    const bool free_sersicindex = true;

    Model::Sources sources{};
    for(size_t i=0; i < 2; ++i)
    {
        std::vector<std::shared_ptr<g2f::Component>> comps;
        for(size_t c=0; c < 2; ++c)
        {
            auto centroids = std::make_shared<g2f::CentroidParameters>(c + i + 1, c + i + 1.5);
            auto integrals = make_integrals(channels, c + 1);
            auto integralmodel = std::make_shared<g2f::LinearIntegralModel>(&integrals);
            std::shared_ptr<g2f::Component> comp;
            if(i == 0)
            {
                comp = std::make_shared<g2f::GaussianComponent>(
                    std::make_shared<g2f::GaussianParametricEllipse>(c + 0.5, c + 1.5, 0),
                    centroids,
                    integralmodel
                );
            } else {
                auto sersic_n = std::make_shared<g2f::SersicMixComponentIndexParameter>(0.5 + 3.5*c);
                sersic_n->set_free(free_sersicindex);
                auto reff_x = std::make_shared<g2f::ReffXParameter>(c + 0.5);
                auto reff_y = std::make_shared<g2f::ReffYParameter>(c + 1.5);
                comp = std::make_shared<g2f::SersicMixComponent>(
                    std::make_shared<g2f::SersicParametricEllipse>(reff_x, reff_y),
                    centroids,
                    integralmodel,
                    sersic_n
                );
            }
            comps.emplace_back(comp);
        }
        auto source = std::make_shared<g2f::Source>(comps);
        sources.push_back(source);
    }

    auto model = std::make_shared<Model>(data, psfmodels, sources);
    CHECK(model->str() != "");

    auto params = model->get_parameters_const_new();
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

    const auto & channel = *channels[0];

    CHECK(model->get_n_gaussians(channel) == 10);
    auto gaussians = model->get_gaussians(channel);
    const size_t n_gauss_src = gaussians->size();
    // 2 comps x 1 (gauss) + 2 comps x 4 (sersic)
    CHECK(n_gauss_src == 10);

    const auto psfcomps = model->get_psfmodels().at(0)->get_gaussians();
    const size_t n_gauss_psf = psfcomps->size();
    CHECK(n_gauss_psf == 2);
    for(size_t p = 0; p < n_gauss_psf; ++p) CHECK(psfcomps->at(p).get_integral_value() > 0);

    auto map_extra = std::make_shared<g2f::extra_param_map>();
    auto map_grad = std::make_shared<g2f::grad_param_map>();

    auto factors_extra = std::make_shared<g2f::extra_param_factors>();
    auto factors_grad = std::make_shared<g2f::grad_param_factors>();

    const size_t n_gauss = n_gauss_src*n_gauss_psf;

    map_extra->reserve(n_gauss);
    map_grad->reserve(n_gauss);
    factors_extra->reserve(n_gauss);
    factors_grad->reserve(n_gauss);

    g2f::ParamCRefs params_src_free_const;
    g2f::ParamRefs params_src_free;
    g2f::ParamFilter filter{false, true, true, true, channel};
    g2f::ParameterMap offsets{};

    for(const auto & source : model->get_sources()) {
        source->get_parameters_const(params_src_free_const, &filter);
        source->get_parameters(params_src_free, &filter);
        for(size_t i_psf = 0; i_psf < n_gauss_psf; ++i_psf)
        {
            source->add_grad_param_map(channel, *map_grad, offsets);
            source->add_extra_param_map(channel, *map_extra, *map_grad, offsets);

            source->add_grad_param_factors(channel, *factors_grad);
            source->add_extra_param_factors(channel, *factors_extra);
        }
    }
    auto n_free = params_src_free.size();
    // 2 sources x (2 comps x (1 integral, 2 centroid, 3 ellipse)) = 24
    // + 1 source x (2 comps x 1 sersicindex) = 26
    CHECK(n_free == 24 + 2*free_sersicindex);

    CHECK(offsets.size() == n_free);

    CHECK(map_extra->size() == n_gauss);
    CHECK(map_grad->size() == n_gauss);
    CHECK(factors_extra->size() == n_gauss);
    CHECK(factors_grad->size() == n_gauss);

    verify_model(*model, channels, params_src_free);
}

TEST_CASE("Model PSF") {
    const Channels channels = CHANNELS_NONE;
    auto data = make_data(channels, 11, 13);
    std::vector<std::shared_ptr<g2f::Component>> comps = {};
    auto integrals = make_integrals(CHANNELS_NONE, 7.0, true);
    auto model_total = std::make_shared<g2f::LinearIntegralModel>(&integrals);
    std::shared_ptr<g2f::IntegralModel> last = model_total;

    const bool fixed = false;
    auto cens = std::make_shared<g2f::CentroidParameters>(
        std::make_shared<g2f::CentroidXParameter>(0, nullptr, nullptr, nullptr, fixed),
        std::make_shared<g2f::CentroidYParameter>(0, nullptr, nullptr, nullptr, fixed)
    );

    for(const auto & sizefrac : {std::pair{1.5, 1.0-1e-15}, {2.5, 1.0}})
    {
        const auto is_last = sizefrac.second == 1;
        auto frac = std::make_shared<g2f::ProperFractionParameter>(sizefrac.second);
        // Would set this condition if we wanted the other fractions free
        if(is_last) frac->set_fixed(true);
        g2f::FractionalIntegralModel::Data data {{g2f::Channel::NONE(), frac}};
        const double & size = sizefrac.first;

        last = g2f::FractionalIntegralModel::make(data, last, is_last);

        comps.emplace_back(std::make_shared<g2f::GaussianComponent>(
            std::make_shared<g2f::GaussianParametricEllipse>(
                std::make_shared<g2f::SigmaXParameter>(size, nullptr, nullptr, nullptr, fixed),
                std::make_shared<g2f::SigmaYParameter>(size, nullptr, nullptr, nullptr, fixed),
                std::make_shared<g2f::RhoParameter>(0, nullptr, nullptr, nullptr, fixed)
            ),
            cens,
            last
        ));
    }
    g2f::Components psfcomps = {g2f::GaussianComponent::make_uniq_default_gaussians({0.})};
    Model::PsfModels psfmodels({std::make_shared<g2f::PsfModel>(psfcomps)});
    Model::Sources sources({std::make_shared<g2f::Source>(comps)});
    auto model = std::make_shared<Model>(data, psfmodels, sources);

    g2f::ParamRefs params_free;
    g2f::ParamFilter filter{false, true, true, true, g2f::Channel::NONE()};
    g2f::ParameterMap offsets{};

    model->get_parameters(params_free, &filter);
    params_free = g2f::nonconsecutive_unique(params_free);
    // 2 cens + 1 fluxfrac + 2*(2sigmas + rho) = 9
    CHECK(params_free.size() == 9);

    verify_model(*model, channels, params_free, true, true);
}
