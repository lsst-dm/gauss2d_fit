#include "ellipticalcomponent.h"
#include "gauss2d/fit/param_defs.h"
#include "param_filter.h"
#include "parameters.h"
#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN

#include "doctest.h"

#include <set>

#include "gauss2d/vectorimage.h"

#include "channel.h"
#include "data.h"
#include "fractionalintegralmodel.h"
#include "gaussiancomponent.h"
#include "integralmodel.h"
#include "linearintegralmodel.h"
#include "model.h"
#include "observation.h"

typedef gauss2d::VectorImage<double> Image;
typedef gauss2d::VectorImage<size_t> Indices;
typedef gauss2d::VectorImage<bool> Mask;
typedef g2f::Observation<double, Image, Mask> Observation;
typedef g2f::Data<double, Image, Mask> Data;
typedef g2f::Model<double, Image, Indices, Mask> Model;

namespace g2f = gauss2d::fit;

const std::vector<std::shared_ptr<const g2f::Channel>> CHANNELS_NONE = {g2f::Channel::NONE_PTR()};

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

TEST_CASE("Model") {
    const std::vector<std::shared_ptr<const g2f::Channel>> channels = {
        // TODO: Figure out how this works - auto const conversion?
        g2f::Channel::make("r"),
        g2f::Channel::make("g"),
        g2f::Channel::make("b")
    };
    const size_t N_X = 5, N_Y = 5;
    std::vector<std::shared_ptr<const Observation>> observations;
    for(const auto & channel_ptr : channels)
    {
        observations.emplace_back(std::make_shared<Observation>(
            std::make_unique<Image>(N_Y, N_X),
            std::make_unique<Image>(N_Y, N_X),
            std::make_unique<Mask>(N_Y, N_X),
            *channel_ptr
        ));
    }
    auto data = std::make_shared<Data>(observations);
    
    Model::PsfModels psfmodels{};
    for(size_t i = 0; i < observations.size(); ++i) {
        g2f::PsfModel::Components comps;
        auto integrals = make_integrals(CHANNELS_NONE, 1.0, true);
        auto model_total =  std::make_shared<g2f::LinearIntegralModel>(&integrals);
        std::shared_ptr<g2f::IntegralModel> last = model_total;
        for(const auto & sizefrac : {std::pair{1, 1/3}, {2, 1}})
        {
            auto frac = std::make_shared<g2f::ProperFractionParameter>(sizefrac.second);
            if(sizefrac.second == 1) frac->set_fixed(true);
            g2f::FractionalIntegralModel::Data data {{g2f::Channel::NONE(), frac}};
            auto model_frac = g2f::FractionalIntegralModel::make(data, last);

            auto integrals = make_integrals(CHANNELS_NONE, sizefrac.second, true);
            comps.emplace_back(std::make_unique<g2f::GaussianComponent>(
                std::make_shared<g2f::GaussianParametricEllipse>(sizefrac.first, sizefrac.first, 0),
                nullptr,
                model_frac
            ));
            last = model_frac;
        }
        psfmodels.push_back(std::make_shared<g2f::PsfModel>(comps));
    }

    g2f::ParamCRefs params_data{};
    data->get_parameters_const(params_data);
    CHECK(params_data.size() == 0);

    Model::Sources sources{};
    for(size_t i=0; i < 2; ++i)
    {
        std::vector<std::shared_ptr<g2f::Component>> comps;
        for(size_t c=0; c < 2; ++c)
        {
            auto integrals = make_integrals(channels, c + 1);
            comps.emplace_back(std::make_unique<g2f::GaussianComponent>(
                std::make_shared<g2f::GaussianParametricEllipse>(c + 1., c + 1., 0),
                nullptr,
                std::make_shared<g2f::LinearIntegralModel>(&integrals)
            ));
        }
        auto source = std::make_shared<g2f::Source>(comps);
        sources.push_back(source);
    }

    auto model = std::make_shared<Model>(data, psfmodels, sources);
    CHECK(model->str() != "");

    auto params = model->get_parameters_const_new();
    // 2 sources x (2 comps x (3 integral, 2 centroid, 3 ellipse)) = 32
    const size_t n_params_src = 32;
    // PSF: 3 observations x (2 comp x (1 integral, n frac, 2 centroid, 3 ellipse)) = 45
    // (each comp after first has an extra frac per channel)
    const size_t n_params_psf = 45;
    CHECK(params.size() == (n_params_src + n_params_psf));

    // PSF: 3 observations x (1 integral + 2 comp x (1 frac, 2 centroid, 3 ellipse)) = 39
    std::set<g2f::ParamBaseCRef> paramset(params.cbegin(), params.cend());
    const size_t n_params_psf_uniq = 39;
    CHECK(paramset.size() == n_params_src + n_params_psf_uniq);

    const auto & channel = *channels[0];

    auto gaussians = model->get_gaussians(channel);
    const size_t n_gauss_src = gaussians->size();
    CHECK(n_gauss_src == 4);

    const auto psfcomps = model->get_psfmodels().at(0)->get_gaussians();
    const size_t n_gauss_psf = psfcomps->size();
    CHECK(n_gauss_psf == 2);

    auto map_extra = std::make_shared<g2f::extra_param_map>();
    auto map_grad = std::make_shared<g2f::grad_param_map>();

    auto factors_extra = std::make_shared<g2f::extra_param_factors>();
    auto factors_grad = std::make_shared<g2f::grad_param_factors>();

    const size_t n_gauss = n_gauss_src*n_gauss_psf;

    map_extra->reserve(n_gauss);
    map_grad->reserve(n_gauss);
    factors_extra->reserve(n_gauss);
    factors_grad->reserve(n_gauss);

    g2f::ParamCRefs params_src_free;
    g2f::ParamFilter filter{false, true, true, true, channel};
    g2f::ParameterMap offsets{};

    for(const auto & source : model->get_sources()) {
        source->get_parameters_const(params_src_free, &filter);
        for(size_t i_psf = 0; i_psf < n_gauss_psf; ++i_psf)
        {
            source->add_grad_param_map(channel, *map_grad, offsets);
            source->add_extra_param_map(channel, *map_extra, offsets);

            source->add_grad_param_factors(channel, *factors_grad);
            source->add_extra_param_factors(channel, *factors_extra);
        }
    }
    auto n_free = params_src_free.size();
    // 2 sources x (2 comps x (1 integral, 2 centroid, 3 ellipse)) = 24
    CHECK(n_free == 24);

    CHECK(offsets.size() == n_free);

    CHECK(map_extra->size() == n_gauss);
    CHECK(map_grad->size() == n_gauss);
    CHECK(factors_extra->size() == n_gauss);
    CHECK(factors_grad->size() == n_gauss);

    for(size_t i = 0; i < 2; ++i)
    {
        for(unsigned short do_gradients=false; do_gradients <= true; do_gradients++)
        {            
            model->setup_evaluators(do_gradients);

            std::vector<double> result = model->evaluate();

            CHECK(result[0] == 0);
            CHECK(result[1] == 0);
            CHECK(result[2] == 0);

            auto outputs = model->get_outputs();
            if(!do_gradients) {
                CHECK(outputs.size() == channels.size());
                CHECK(outputs[0]->get_value(0, 0) == outputs[1]->get_value(0, 0));
                CHECK(outputs[0]->get_value(0, 0) == outputs[2]->get_value(0, 0));
            }
        }
    }
}