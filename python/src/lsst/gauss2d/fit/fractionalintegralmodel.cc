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

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <memory>

#include "lsst/gauss2d/fit/data.h"
#include "lsst/gauss2d/fit/integralmodel.h"
#include "lsst/gauss2d/fit/fractionalintegralmodel.h"
#include "lsst/gauss2d/fit/parameters.h"
#include "pybind11.h"

namespace py = pybind11;
using namespace pybind11::literals;

namespace g2f = lsst::gauss2d::fit;

void bind_fractionalintegralmodel(py::module &m) {
    auto _p = py::class_<g2f::FractionalIntegralModel, std::shared_ptr<g2f::FractionalIntegralModel>,
                         g2f::IntegralModel>(m, "FractionalIntegralModel")
                      .def(py::init(&g2f::FractionalIntegralModel::make), "data"_a, "model"_a,
                           "is_final"_a = false)
                      .def_property_readonly("channels", &g2f::FractionalIntegralModel::get_channels)
                      .def("integral", &g2f::FractionalIntegralModel::get_integral, "channel"_a)
                      .def("integral_remainder", &g2f::FractionalIntegralModel::get_integral_remainder,
                           "channel"_a)
                      .def("parameters", &g2f::FractionalIntegralModel::get_parameters,
                           "parameters"_a = g2f::ParamRefs(), "paramfilter"_a = nullptr)
                      .def_property_readonly("parent_model", &g2f::FractionalIntegralModel::get_parent_model)
                      .def_static("model", &g2f::FractionalIntegralModel::find_model)
                      .def("__getitem__", [](const g2f::FractionalIntegralModel &self,
                                             const g2f::Channel &c) { return self.at(c); })
                      .def("__len__", &g2f::FractionalIntegralModel::size)
                      .def("__repr__",
                           [](const g2f::FractionalIntegralModel &self) { return self.repr(true); })
                      .def("__str__", &g2f::FractionalIntegralModel::str);
}
