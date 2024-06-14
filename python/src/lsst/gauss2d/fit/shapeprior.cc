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

#include "lsst/gauss2d/fit/parametricellipse.h"
#include "lsst/gauss2d/fit/parametricgaussian1d.h"
#include "lsst/gauss2d/fit/prior.h"
#include "lsst/gauss2d/fit/shapeprior.h"
#include "pybind11.h"

namespace py = pybind11;
using namespace pybind11::literals;

namespace g2f = lsst::gauss2d::fit;

void bind_shapeprior(py::module &m) {
    auto _o = py::class_<g2f::ShapePriorOptions, std::shared_ptr<g2f::ShapePriorOptions>, lsst::gauss2d::Object>(
                      m, "ShapePriorOptions")
                      .def(py::init<double, double, double>(),
                           "delta_jacobian"_a = g2f::ShapePriorOptions::delta_jacobian_default,
                           "size_maj_floor"_a = g2f::ShapePriorOptions::size_maj_floor_default,
                           "axrat_floor"_a = g2f::ShapePriorOptions::axrat_floor_default)
                      .def("check_delta_jacobian", &g2f::ShapePriorOptions::check_delta_jacobian)
                      .def("check_size_maj_floor", &g2f::ShapePriorOptions::check_size_maj_floor)
                      .def("check_axrat_floor", &g2f::ShapePriorOptions::check_axrat_floor)
                      .def_property("delta_jacobian", &g2f::ShapePriorOptions::get_delta_jacobian,
                                    &g2f::ShapePriorOptions::set_delta_jacobian)
                      .def_property("size_maj_floor", &g2f::ShapePriorOptions::get_size_maj_floor,
                                    &g2f::ShapePriorOptions::set_size_maj_floor)
                      .def_property("size_maj_floor", &g2f::ShapePriorOptions::get_size_maj_floor,
                                    &g2f::ShapePriorOptions::set_size_maj_floor)
                      .def("__repr__", [](const g2f::ShapePriorOptions &self) { return self.repr(true); })
                      .def("__str__", &g2f::ShapePriorOptions::str);

    auto _s = py::class_<g2f::ShapePrior, std::shared_ptr<g2f::ShapePrior>, g2f::Prior>(m, "ShapePrior")
                      .def(py::init<std::shared_ptr<const g2f::ParametricEllipse>,
                                    std::shared_ptr<g2f::ParametricGaussian1D>,
                                    std::shared_ptr<g2f::ParametricGaussian1D>,
                                    std::shared_ptr<g2f::ShapePriorOptions> >(),
                           "ellipse"_a, "prior_size"_a = nullptr, "prior_axrat"_a = nullptr,
                           "options"_a = nullptr)
                      .def("evaluate", &g2f::ShapePrior::evaluate, "calc_jacobians"_a = false,
                           "normalize_loglike"_a = false)
                      .def_property_readonly("loglike_const_terms", &g2f::ShapePrior::get_loglike_const_terms)
                      .def_property("prior_size", &g2f::ShapePrior::get_prior_size,
                                    &g2f::ShapePrior::set_prior_size)
                      .def_property("prior_axrat", &g2f::ShapePrior::get_prior_axrat,
                                    &g2f::ShapePrior::set_prior_axrat)
                      .def("__len__", &g2f::ShapePrior::size)
                      .def("__repr__", [](const g2f::ShapePrior &self) { return self.repr(true); })
                      .def("__str__", &g2f::ShapePrior::str);
}
