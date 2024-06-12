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
#include "lsst/gauss2d/fit/sersicparametricellipse.h"
#include "pybind11.h"

namespace py = pybind11;
using namespace pybind11::literals;

namespace g2f = lsst::gauss2d::fit;

void bind_sersicparametricellipse(py::module &m) {
    auto _e = py::class_<g2f::SersicParametricEllipse, std::shared_ptr<g2f::SersicParametricEllipse>,
                         g2f::ParametricEllipse>(m, "SersicParametricEllipse")
                      .def(py::init<double, double, double>(), "size_x"_a = 0, "size_y"_a = 0, "rho"_a = 0)
                      .def(py::init<std::shared_ptr<g2f::ReffXParameterD>,
                                    std::shared_ptr<g2f::ReffYParameterD>,
                                    std::shared_ptr<g2f::RhoParameterD> >(),
                           "size_x"_a = nullptr, "size_y"_a = nullptr, "rho"_a = nullptr)
                      .def_property("rho", &g2f::SersicParametricEllipse::get_rho,
                                    &g2f::SersicParametricEllipse::set_rho)
                      .def_property("size_x", &g2f::SersicParametricEllipse::get_size_x,
                                    &g2f::SersicParametricEllipse::set_size_x)
                      .def_property("size_y", &g2f::SersicParametricEllipse::get_size_y,
                                    &g2f::SersicParametricEllipse::set_size_y)
                      .def_property("xyr", &g2f::SersicParametricEllipse::get_xyr,
                                    &g2f::SersicParametricEllipse::set_xyr)
                      .def_property_readonly("rho_param", &g2f::SersicParametricEllipse::get_rho_param)
                      .def_property_readonly("size_x_param", &g2f::SersicParametricEllipse::get_size_x_param)
                      .def_property_readonly("size_y_param", &g2f::SersicParametricEllipse::get_size_y_param)
                      .def_property_readonly("rho_param_ptr",
                                             &g2f::SersicParametricEllipse::get_rho_param_ptr)
                      .def_property_readonly("size_x_param_ptr",
                                             &g2f::SersicParametricEllipse::get_size_x_param_ptr)
                      .def_property_readonly("size_y_param_ptr",
                                             &g2f::SersicParametricEllipse::get_size_y_param_ptr)
                      .def("__repr__",
                           [](const g2f::SersicParametricEllipse &self) { return self.repr(true); })
                      .def("__str__", &g2f::SersicParametricEllipse::str);
}
