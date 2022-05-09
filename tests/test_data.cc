#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN

#include "doctest.h"

#include <memory>
#include <stdexcept>

#include "gauss2d/vectorimage.h"

#include "data.h"
#include "observation.h"

namespace g2f = gauss2d::fit;

typedef gauss2d::VectorImage<double> Image;
typedef gauss2d::VectorImage<bool> Mask;
typedef g2f::Observation<double, Image, Mask> Observation;
typedef g2f::Data<double, Image, Mask> Data;

TEST_CASE("Data") {
    const size_t N_X = 5, N_Y = 5;
    std::vector<std::shared_ptr<const Observation>> observations = {
        std::make_shared<Observation>(
            std::make_unique<Image>(N_Y, N_X),
            std::make_unique<Image>(N_Y, N_X),
            std::make_unique<Mask>(N_Y, N_X)
        ),
        std::make_shared<Observation>(
            std::make_unique<Image>(N_Y, N_X + 1),
            std::make_unique<Image>(N_Y, N_X + 1),
            std::make_unique<Mask>(N_Y, N_X + 1)
        )
    };
    const Observation & observation = *(observations.at(0));

    auto data = std::make_shared<Data>(observations);
    g2f::ParamCRefs params{};
    // no parameters yet (until background model added)
    CHECK(data->get_parameters_const(params).size() == 0);
    const Observation & observation_data = data->at(0).get();
    CHECK(&observation == &observation_data);
    CHECK(observation == observation_data);
    CHECK(data->str() != "");
}