#include "ellipticalcomponent.h"
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
        auto integrals = make_integrals(CHANNELS_NONE, 1.0, true);
        g2f::PsfModel::Components components;
        components.emplace_back(std::make_unique<g2f::GaussianComponent>(
            nullptr,
            std::make_shared<g2f::EllipseParameters>(1., 1., 0),
            std::make_shared<g2f::LinearIntegralModel>(&integrals)
        ));
        psfmodels.push_back(std::make_shared<g2f::PsfModel>(components));
    }

    g2f::ParamCRefs params_data{};
    data->get_parameters_const(params_data);
    CHECK(params_data.size() == 0);

    Model::Sources sources{};
    for(size_t i=0; i < 2; ++i)
    {
        std::vector<std::shared_ptr<g2f::Component>> comps;
        std::shared_ptr<g2f::IntegralModel> last;
        for(size_t c=0; c < 2; ++c)
        {
            if(c == 0) {
                g2f::LinearIntegralModel::Data integrals{};
                g2f::FractionalIntegralModel::Data fracs{};
                for(const auto & channel: channels) {
                    integrals[*channel] = std::make_shared<g2f::IntegralParameter>(c + 1.);
                    fracs[*channel] = std::make_shared<g2f::ProperFractionParameter>(0.5);
                }
                auto integralmodel = std::make_shared<const g2f::LinearIntegralModel>(&integrals);
                last = g2f::FractionalIntegralModel::make(fracs, integralmodel);
            } else {
                g2f::FractionalIntegralModel::Data fracs{};
                for(const auto & channel: channels) {
                    auto frac = std::make_shared<g2f::ProperFractionParameter>(1.0);
                    frac->set_fixed(true);
                    fracs[*channel] = frac;
                }
                last = g2f::FractionalIntegralModel::make(fracs, last);
            }

            auto comp = std::make_unique<g2f::GaussianComponent>(
                nullptr,
                std::make_shared<g2f::EllipseParameters>(c + 1., c + 1., 0),
                last
            );
            comps.push_back(std::move(comp));
        }
        auto source = std::make_shared<g2f::Source>(comps);
        sources.push_back(source);
    }

    auto model = std::make_shared<Model>(data, psfmodels, sources);
    CHECK(model->str() != "");

    auto params = model->get_parameters_const_new();
    // 2 comps x (3 integral, 3 frac, 2 centroid, 3 ellipse) = 22
    // 2 comps x (3 integral, 6 frac, 2 centroid, 3 ellipse) = 28
    // (each comp after first has an extra frac per channel)
    // PSF: 3 observations x (1 comp x (1 integral, 2 centroid, 3 ellipse)) = 18
    CHECK(params.size() == 68);

    // 2 comps x (3 integral, 3 frac, 2 centroid, 3 ellipse) = 22
    // 2 comps x (3 frac, 2 centroid, 3 ellipse) = 16
    // PSF: 3 observations x (1 comp x (1 integral, 2 centroid, 3 ellipse)) = 18
    std::set<g2f::ParamBaseCRef> paramset(params.cbegin(), params.cend());
    CHECK(paramset.size() == 56);

    auto gaussians = model->get_gaussians(*channels[0]);
    CHECK(gaussians->size() == 4);

    for(size_t i = 0; i < 3; ++i)
    {
        model->setup_evaluators();

        std::vector<double> result = model->evaluate();

        CHECK(result[0] == 0);
        CHECK(result[1] == 0);
        CHECK(result[2] == 0);

        auto outputs = model->get_outputs();
        CHECK(outputs.size() == channels.size());
        CHECK(outputs[0]->get_value(0, 0) == outputs[1]->get_value(0, 0));
        CHECK(outputs[0]->get_value(0, 0) == outputs[2]->get_value(0, 0));
    }
}