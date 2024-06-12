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

#include <memory>

#include <pybind11/attr.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "pybind11.h"

#include "lsst/gauss2d/fit/component.h"
#include "lsst/gauss2d/fit/gaussiancomponent.h"
#include "lsst/gauss2d/fit/gaussianparametricellipse.h"

namespace py = pybind11;
using namespace pybind11::literals;

namespace g2f = lsst::gauss2d::fit;

void bind_gaussiancomponent(py::module &m) {
    auto _e = py::class_<g2f::GaussianComponent, std::shared_ptr<g2f::GaussianComponent>,
                         g2f::EllipticalComponent>(m, "GaussianComponent")
                      .def(py::init<std::shared_ptr<g2f::GaussianParametricEllipse>,
                                    std::shared_ptr<g2f::CentroidParameters>,
                                    std::shared_ptr<g2f::IntegralModel>>(),
                           "ellipse"_a = nullptr, "centroid"_a = nullptr, "integral"_a = nullptr)
                      .def("parameters", &g2f::GaussianComponent::get_parameters,
                           "parameters"_a = g2f::ParamRefs(), "paramfilter"_a = nullptr)
                      .def("gaussians",
                           [](const g2f::GaussianComponent &g, const g2f::Channel &c) {
                               return std::shared_ptr<const lsst::gauss2d::Gaussians>(g.get_gaussians(c));
                           })
                      .def_static("make_uniq_default_gaussians",
                                  &g2f::GaussianComponent::make_uniq_default_gaussians, "sizes"_a, "fixed"_a)
                      .def("__repr__", [](const g2f::GaussianComponent &self) { return self.repr(true); })
                      .def("__str__", &g2f::GaussianComponent::str);
}
