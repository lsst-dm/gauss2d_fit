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

#include "gauss2d/fit/param_filter.h"
#include "pybind11.h"

namespace py = pybind11;
using namespace pybind11::literals;

namespace g2f = gauss2d::fit;

void bind_param_filter(py::module &m)
{
    auto _p = py::class_<g2f::ParamFilter, std::shared_ptr<g2f::ParamFilter>>(m, "ParamFilter")
        .def(py::init<bool, bool, bool, bool>(), "fixed"_a, "free"_a, "linear"_a, "nonlinear"_a)
        .def_readwrite("fixed", &g2f::ParamFilter::fixed)
        .def_readwrite("free", &g2f::ParamFilter::free)
        .def_readwrite("linear", &g2f::ParamFilter::linear)
        .def_readwrite("nonlinear", &g2f::ParamFilter::nonlinear)
    ;
}
