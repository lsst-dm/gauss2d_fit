/*
 * This file is part of gauss2dfit.
 *
 * Developed for the LSST Data Management System.
 * This product includes software developed by the LSST Project
 * (https://www.lsst.org).
 * See the COPYRIGHT file at the top-level directory of this distribution
 * for details of code ownership.
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */

#ifdef GAUSS2D_FIT_HAS_GSL

#ifndef GAUSS2D_FIT_GSLIMAGE_H
#define GAUSS2D_FIT_GSLIMAGE_H

#include "gsl.h"

#include "gauss2d/object.h"

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_matrix_float.h>
#include <gsl/gsl_matrix_ulong.h>

#include <memory>
#include <stdexcept>

#include "gauss2d/evaluate.h"
#include "gauss2d/image.h"

namespace gauss2d::fit {

#pragma GCC visibility push(hidden)

template <typename T>
struct GSLMatrix {};

template <typename T, typename matrix_class, matrix_class *alloc_func(const size_t, const size_t),
          void free_func(matrix_class *)>
struct GSLMatrixImpl {
    using matrix_type = matrix_class;
    matrix_type *matrix = nullptr;

    inline matrix_class *alloc(const size_t n1, const size_t n2) { return alloc_func(n1, n2); }

    inline void free() { free_func(matrix); }

    inline T get_unchecked(const size_t i, const size_t j) const {
        // It's easier to just copy this then fiddle with defines to
        // disable GSL's range checking while calling gsl_matrix_get
        return matrix->data[i * matrix->tda + j];
    }

    inline T *ptr(const size_t i, const size_t j) { return (T *)(matrix->data + (i * matrix->tda + j)); }
};

// TODO: Shouldn't this be possible with using or typedef or something instead of inheritance?
template <>
struct GSLMatrix<double> : GSLMatrixImpl<double, gsl_matrix, gsl_matrix_alloc, gsl_matrix_free> {};

template <>
struct GSLMatrix<float>
        : GSLMatrixImpl<float, gsl_matrix_float, gsl_matrix_float_alloc, gsl_matrix_float_free> {};

template <>
struct GSLMatrix<size_t>
        : GSLMatrixImpl<unsigned long, gsl_matrix_ulong, gsl_matrix_ulong_alloc, gsl_matrix_ulong_free> {};

template <typename t>
class GSLImage : public gauss2d::Image<t, GSLImage<t>> {
public:
    using data = GSLMatrix<t>;
    typename data::matrix_type matrix_type;

private:
    const size_t _n_rows;
    const size_t _n_cols;

    data _data;

public:
    inline t &_get_value_unchecked_impl(size_t row, size_t col) { return *(_data.ptr(row, col)); }

    size_t get_n_cols_impl() const { return _n_cols; };
    size_t get_n_rows_impl() const { return _n_rows; };

    const inline t get_value_unchecked_impl(size_t row, size_t col) const {
        return _data.get_unchecked(row, col);
    }

    /*
    void set_value_unchecked(size_t row, size_t col, t value) {
        this->_data_ref(row, col) = value;
    }*/

    GSLImage(size_t n_rows, size_t n_cols, const t *value_init = Image<t, GSLImage<t>>::_value_default_ptr(),
             const std::shared_ptr<const gauss2d::CoordinateSystem> coordsys = nullptr)
            : gauss2d::Image<t, GSLImage<t>>(coordsys), _n_rows(n_rows), _n_cols(n_cols) {
        _data.matrix = _data.alloc(_n_rows, _n_cols);
        if (value_init != nullptr) {
            this->fill(*value_init);
        }
    }

    ~GSLImage() {
        if (_data.matrix != nullptr) {
            _data.free();
        }
    }
};
#pragma GCC visibility pop

/*
template<typename T>
void declare_image(py::module &m, std::string typestr) {
    using Class = GSLImage<T>;
    std::string pyclass_name = std::string("Image") + typestr;
    py::class_<Class, std::shared_ptr<Class>>(m, pyclass_name.c_str())
    .def(
        py::init<size_t, size_t, const std::shared_ptr<const gauss2d::CoordinateSystem>>(),
        "n_rows"_a, "n_cols"_a, "coordsys"_a = gauss2d::COORDS_DEFAULT
    )
    .def(
        py::init<py::array_t<T>, const std::shared_ptr<const gauss2d::CoordinateSystem>>(),
        "data"_a, "coordsys"_a = gauss2d::COORDS_DEFAULT
    )
    .def_property_readonly("n_rows", &Class::get_n_rows)
    .def_property_readonly("n_cols", &Class::get_n_cols)
    .def_property_readonly("data", &Class::get_data)
    .def("fill", &Class::fill)
    .def("get_value", &Class::get_value)
    .def("set_value", &Class::set_value)
    .def("get_value_unchecked", &Class::get_value_unchecked)
    .def("set_value_unchecked", &Class::set_value_unchecked)
    .def_property_readonly("size", &Class::size)
    .def_property_readonly("shape", &Class::shape)
    .def("__repr__", [](const Class & self) { return self.repr(true); })
    .def("__str__", &Class::str)
    ;
}

template<typename T>
void declare_image_array(py::module &m, std::string typestr) {
    using Class = gauss2d::ImageArray<T, GSLImage<T>>;
    std::string pyclass_name = std::string("ImageArray") + typestr;
    py::class_<Class, std::shared_ptr<Class>>(m, pyclass_name.c_str())
    .def(py::init<const typename Class::Data*>(), "data"_a)
    .def("at", &Class::at, py::return_value_policy::reference)
    .def_property_readonly("size", &Class::size)
    .def("__repr__", [](const Class & self) { return self.repr(true); })
    .def("__str__", &Class::str)
    ;
}

template<typename t>
void declare_evaluator(py::module &m, std::string typestr) {
    using Class = gauss2d::GaussianEvaluator<t, GSLImage<t>, GSLImage<gauss2d::idx_type>>;
    std::string pyclass_name = std::string("GaussianEvaluator") + typestr;
    py::class_<Class, std::shared_ptr<Class>>(m, pyclass_name.c_str())
    .def(
        py::init<
            const std::shared_ptr<const gauss2d::ConvolvedGaussians>,
            const std::shared_ptr<const gauss2d::CoordinateSystem>,
            const std::shared_ptr<const GSLImage<t>>,
            const std::shared_ptr<const GSLImage<t>>,
            const std::shared_ptr<GSLImage<t>>,
            const std::shared_ptr<GSLImage<t>>,
            const std::shared_ptr<gauss2d::ImageArray<t, GSLImage<t>>>,
            const std::shared_ptr<const GSLImage<gauss2d::idx_type>>,
            const std::shared_ptr<const GSLImage<t>>,
            const std::shared_ptr<const GSLImage<gauss2d::idx_type>>,
            const std::shared_ptr<const GSLImage<t>>,
            const std::shared_ptr<const GSLImage<t>>
        >(),
        "gaussians"_a,
        "coordsys"_a = nullptr,
        "data"_a = nullptr,
        "sigma_inv"_a = nullptr,
        "output"_a = nullptr,
        "residual"_a = nullptr,
        "grads"_a = nullptr,
        "grad_param_map"_a = nullptr,
        "grad_param_factor"_a = nullptr,
        "extra_param_map"_a = nullptr,
        "extra_param_factor"_a = nullptr,
        "background"_a = nullptr
    )
    .def("loglike_pixel", &Class::loglike_pixel, "to_add"_a = false)
    .def_property_readonly("n_cols", &Class::get_n_cols)
    .def_property_readonly("n_rows", &Class::get_n_rows)
    .def_property_readonly("n_size", &Class::get_size)
    ;
}

template<typename t, class Data, class Indices>
void declare_maker(py::module &m, std::string typestr) {
    m.def(
        ("make_gaussians_pixel_" + typestr).c_str(), gauss2d::make_gaussians_pixel<t, Data, Indices>,
        "Evaluate a 2D Gaussian at the centers of pixels on a rectangular grid using the standard bivariate"
        "Gaussian PDF.",
        "gaussians"_a, "output"_a=nullptr, "n_rows"_a=0, "n_cols"_a=0, "coordsys"_a=nullptr
    );
}
*/

}  // namespace gauss2d::fit

#endif  // GAUSS2D_FIT_GSLIMAGE_H
#endif  // GAUSS2D_FIT_HAS_GSL
