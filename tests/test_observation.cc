#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN

#include "doctest.h"

#include <memory>
#include <stdexcept>

#include "lsst/gauss2d/image.h"
#include "lsst/gauss2d/vectorimage.h"

#include "lsst/gauss2d/fit/observation.h"

namespace g2d = lsst::gauss2d;
namespace g2f = lsst::gauss2d::fit;

typedef g2d::VectorImage<double> Image;
typedef g2d::VectorImage<bool> Mask;
typedef g2f::Observation<double, Image, Mask> Observation;

TEST_CASE("Observation") {
    const size_t N_X = 5, N_Y = 5;
    CHECK_THROWS_AS(std::make_shared<Observation>(std::make_unique<Image>(N_Y, N_X + 1),
                                                  std::make_unique<Image>(N_Y, N_X),
                                                  std::make_unique<Mask>(N_Y, N_X)),
                    std::invalid_argument);
    CHECK_THROWS_AS(std::make_shared<Observation>(std::make_unique<Image>(N_Y, N_X),
                                                  std::make_unique<Image>(N_Y, N_X + 1),
                                                  std::make_unique<Mask>(N_Y, N_X)),
                    std::invalid_argument);
    CHECK_THROWS_AS(std::make_shared<Observation>(std::make_unique<Image>(N_Y, N_X),
                                                  std::make_unique<Image>(N_Y, N_X),
                                                  std::make_unique<Mask>(N_Y, N_X + 1)),
                    std::invalid_argument);
    auto observation = std::make_shared<Observation>(std::make_unique<Image>(N_Y, N_X),
                                                     std::make_unique<Image>(N_Y, N_X),
                                                     std::make_unique<Mask>(N_Y, N_X));
    CHECK(lsst::gauss2d::images_compatible<double, Image, double, Image>(observation->get_image(),
                                                                         observation->get_sigma_inverse()));
    CHECK(lsst::gauss2d::images_compatible<double, Image, bool, Mask>(observation->get_image(),
                                                                      observation->get_mask_inverse()));
    CHECK_NE(observation->str(), "");
    g2f::ParamCRefs params{};
    // no parameters yet
    CHECK_EQ(0, observation->get_parameters_const(params).size());
}