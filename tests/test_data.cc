#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN

#include "doctest.h"

#include <memory>
#include <stdexcept>

#include "lsst/gauss2d/vectorimage.h"

#include "lsst/gauss2d/fit/data.h"
#include "lsst/gauss2d/fit/observation.h"

namespace g2f = lsst::gauss2d::fit;

typedef lsst::gauss2d::VectorImage<double> ImageD;
typedef lsst::gauss2d::VectorImage<bool> Mask;
typedef g2f::Observation<double, ImageD, Mask> ObservationD;
typedef g2f::Data<double, ImageD, Mask> DataD;

TEST_CASE("Data") {
    const size_t N_X = 5, N_Y = 5;
    std::vector<std::shared_ptr<const ObservationD>> observations
            = {std::make_shared<ObservationD>(std::make_unique<ImageD>(N_Y, N_X),
                                              std::make_unique<ImageD>(N_Y, N_X),
                                              std::make_unique<Mask>(N_Y, N_X)),
               std::make_shared<ObservationD>(std::make_unique<ImageD>(N_Y, N_X + 1),
                                              std::make_unique<ImageD>(N_Y, N_X + 1),
                                              std::make_unique<Mask>(N_Y, N_X + 1))};
    const ObservationD& observation = *(observations.at(0));

    auto data = std::make_shared<DataD>(observations);
    g2f::ParamCRefs params{};
    // no parameters yet (until background model added)
    CHECK_EQ(data->get_parameters_const(params).size(), 0);
    const ObservationD& observation_data = data->at(0).get();
    CHECK_EQ(&observation, &observation_data);
    CHECK_EQ(observation, observation_data);
    CHECK_NE(data->str(), "");
}