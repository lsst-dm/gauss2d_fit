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

#include "lsst/gauss2d/fit/gaussianprior.h"
#include "lsst/gauss2d/fit/param_defs.h"
#include "lsst/gauss2d/fit/prior.h"
#include "pybind11.h"

namespace py = pybind11;
using namespace pybind11::literals;

namespace g2f = lsst::gauss2d::fit;

void bind_gaussianprior(py::module &m) {
    auto _e = py::class_<g2f::GaussianPrior, std::shared_ptr<g2f::GaussianPrior>, g2f::Prior>(m,
                                                                                              "GaussianPrior")
                      .def(py::init<std::shared_ptr<const g2f::ParamBase>, double, double, bool>(),
                           "param"_a = nullptr, "mean"_a = 0, "stddev"_a = 1., "transformed"_a = false)
                      .def("evaluate", &g2f::GaussianPrior::evaluate, "calc_jacobians"_a = false,
                           "normalize_loglike"_a = false)
                      .def_property_readonly("loglike_const_terms",
                                             &g2f::GaussianPrior::get_loglike_const_terms)
                      .def_property("mean", &g2f::GaussianPrior::get_mean, &g2f::GaussianPrior::set_mean)
                      .def_property_readonly("param", &g2f::GaussianPrior::get_param)
                      .def_property("stddev", &g2f::GaussianPrior::get_stddev,
                                    &g2f::GaussianPrior::set_stddev)
                      .def_property("transformed", &g2f::GaussianPrior::get_transformed,
                                    &g2f::GaussianPrior::set_transformed)
                      .def("__len__", &g2f::GaussianPrior::size)
                      .def("__repr__", [](const g2f::GaussianPrior &self) { return self.repr(true); })
                      .def("__str__", &g2f::GaussianPrior::str);
}
