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

#include <pybind11/attr.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <memory>

#include "lsst/gauss2d/fit/parameters.h"
#include "lsst/gauss2d/fit/parametricgaussian1d.h"

#include "pybind11.h"

namespace py = pybind11;
using namespace pybind11::literals;

namespace g2f = lsst::gauss2d::fit;
void bind_parametricgaussian1d(py::module &m) {
    auto _e = py::class_<g2f::ParametricGaussian1D, std::shared_ptr<g2f::ParametricGaussian1D>,
                         lsst::gauss2d::Object>(m, "ParametricGaussian1D")
                      .def(py::init<std::shared_ptr<g2f::MeanParameterD>,
                                    std::shared_ptr<g2f::StdDevParameterD>>(),
                           "mean"_a = nullptr, "stddev"_a = nullptr)
                      .def_property("mean", &g2f::ParametricGaussian1D::get_mean,
                                    &g2f::ParametricGaussian1D::set_mean)
                      .def_property_readonly("mean_parameter", &g2f::ParametricGaussian1D::get_mean_parameter)
                      .def_property("stddev", &g2f::ParametricGaussian1D::get_stddev,
                                    &g2f::ParametricGaussian1D::set_stddev)
                      .def_property_readonly("stddev_parameter",
                                             &g2f::ParametricGaussian1D::get_stddev_parameter)
                      .def("__repr__", [](const g2f::ParametricGaussian1D &self) { return self.repr(true); })
                      .def("__str__", &g2f::ParametricGaussian1D::str);
}
