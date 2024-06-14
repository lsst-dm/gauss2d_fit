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
#include <string>

#include "lsst/gauss2d/centroid.h"
#include "lsst/gauss2d/fit/centroidparameters.h"
#include "pybind11.h"

// #include "utils.h"

namespace py = pybind11;
using namespace pybind11::literals;

namespace g2f = lsst::gauss2d::fit;

void bind_centroidparameters(py::module &m) {
    auto _c = py::class_<g2f::CentroidParameters, std::shared_ptr<g2f::CentroidParameters>,
                         lsst::gauss2d::CentroidData>(m, "CentroidParameters")
                      .def(py::init<double, double>(), "x"_a = 0, "y"_a = 0)
                      .def(py::init<std::shared_ptr<g2f::CentroidXParameterD>,
                                    std::shared_ptr<g2f::CentroidYParameterD>>(),
                           "x"_a = nullptr, "y"_a = nullptr)
                      .def_property("x", &g2f::CentroidParameters::get_x, &g2f::CentroidParameters::set_x)
                      .def_property("y", &g2f::CentroidParameters::get_y, &g2f::CentroidParameters::set_y)
                      .def_property("xy", &g2f::CentroidParameters::get_xy, &g2f::CentroidParameters::set_xy)
                      .def_property_readonly("x_param", &g2f::CentroidParameters::get_x_param)
                      .def_property_readonly("y_param", &g2f::CentroidParameters::get_y_param)
                      .def_property_readonly("x_param_ptr", &g2f::CentroidParameters::get_x_param_ptr)
                      .def_property_readonly("y_param_ptr", &g2f::CentroidParameters::get_y_param_ptr)
                      .def("__repr__", [](const g2f::CentroidParameters &self) { return self.repr(true); })
                      .def("__str__", &g2f::CentroidParameters::str);
}
