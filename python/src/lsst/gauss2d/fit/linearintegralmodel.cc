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
#include "lsst/gauss2d/fit/linearintegralmodel.h"
#include "lsst/gauss2d/fit/parameters.h"
#include "pybind11.h"

namespace py = pybind11;
using namespace pybind11::literals;

namespace g2f = lsst::gauss2d::fit;

void bind_linearintegralmodel(py::module &m) {
    auto _p = py::class_<g2f::LinearIntegralModel, std::shared_ptr<g2f::LinearIntegralModel>,
                         g2f::IntegralModel>(m, "LinearIntegralModel")
                      .def(py::init<const g2f::LinearIntegralModel::Data *>(), "data"_a)
                      .def_property_readonly("channels", &g2f::LinearIntegralModel::get_channels)
                      .def("integral", &g2f::LinearIntegralModel::get_integral, "channel"_a)
                      .def("parameters", &g2f::LinearIntegralModel::get_parameters,
                           "parameters"_a = g2f::ParamRefs(), "paramfilter"_a = nullptr)
                      .def("__getitem__", [](const g2f::LinearIntegralModel &self,
                                             const g2f::Channel &c) { return self.at(c); })
                      .def("__len__", &g2f::LinearIntegralModel::size)
                      .def("__repr__", [](const g2f::LinearIntegralModel &self) { return self.repr(true); })
                      .def("__str__", &g2f::LinearIntegralModel::str);
    /*
        typename Data::iterator begin() noexcept;
        typename Data::const_iterator cbegin() const noexcept;

        typename Data::iterator end() noexcept;
        typename Data::const_iterator cend() const noexcept;
    */
}
