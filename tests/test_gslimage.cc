#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN

#include "doctest.h"

#ifdef GAUSS2D_FIT_HAS_GSL

#include "gsl.h"
#include "gslimage.h"

namespace g2 = gauss2d;

typedef double number;

typedef g2::fit::GSLImage<number> Image;
typedef g2::ImageArray<number, Image> ImageArray;

TEST_CASE("VectorImage") {
    const number value_init = -42.1;
    Image nonzero{1, 1, &value_init};
    CHECK(nonzero.get_value(0, 0) == value_init);

    size_t n_rows = 3, n_cols = 2;
    std::shared_ptr<Image> image = std::make_shared<Image>(n_rows, n_cols);
    Image zeros{n_rows, n_cols};

    auto x = Image(n_rows, n_cols) == zeros;
    CHECK(x);

    CHECK(image->get_coordsys() == g2::COORDS_DEFAULT);
    CHECK(&(image->get_coordsys()) == &(g2::COORDS_DEFAULT));

    CHECK(image->get_n_cols() == n_cols);
    CHECK(image->get_n_rows() == n_rows);
    CHECK(image->size() == n_cols * n_rows);

    image->add_value_unchecked(0, 0, 1);
    CHECK(*image != zeros);
    CHECK(image->get_value_unchecked(0, 0) == 1);
    CHECK_THROWS_AS(image->get_value(n_rows, n_cols), std::out_of_range);
    CHECK_THROWS_AS(image->set_value(n_rows, n_cols, 1), std::out_of_range);
    auto & value = image->_get_value_unchecked(1, 1);
    value = -1;
    CHECK(image->get_value_unchecked(1, 1) == -1);
}

TEST_CASE("VectorImageArray") {
    size_t n_rows = 3, n_cols = 2;
    std::shared_ptr<Image> image = std::make_shared<Image>(n_rows, n_cols);
    ImageArray::Data data = {image, image};
    ImageArray arr = ImageArray(&data);

    for (ImageArray::const_iterator it = arr.cbegin(); it != arr.cend(); ++it) {
        CHECK(*it == image);
    }
    for (const auto& img : arr) {
        CHECK(img == image);
    }
}

#endif