/*
 * This file is part of gauss2dfit.
 *
 * Developed for the LSST Source Management System.
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

#include "gauss2d/fit/parametricmodel.h"
#include "gauss2d/fit/source.h"
#include "pybind11.h"

namespace py = pybind11;
using namespace pybind11::literals;

namespace g2f = gauss2d::fit;

void bind_source(py::module &m)
{
    auto _s = py::class_<g2f::Source,
        std::shared_ptr<g2f::Source>,
        g2f::ParametricModel
    >(m, "Source")
        .def(py::init<g2f::Source::Components &>(), "components"_a=nullptr)
        .def("gaussians", [](const g2f::Source &g, const g2f::Channel & c)
            { return std::shared_ptr<gauss2d::Gaussians>(g.get_gaussians(c)); })
        .def("parameters", &g2f::Source::get_parameters, "parameters"_a=g2f::ParamRefs(), "paramfilter"_a=nullptr)
        .def("__repr__", &g2f::Source::str)
    ;
}
