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

#include "lsst/gauss2d/object.h"
#include "lsst/gauss2d/fit/param_defs.h"
#include "lsst/gauss2d/fit/prior.h"
#include "pybind11.h"

namespace py = pybind11;
using namespace pybind11::literals;

namespace g2f = lsst::gauss2d::fit;

void bind_prior(py::module &m) {
    auto _e = py::class_<g2f::PriorEvaluation, std::shared_ptr<g2f::PriorEvaluation>, lsst::gauss2d::Object>(
                      m, "PriorEvaluation")
                      .def(py::init<double, std::vector<double>,
                                    std::map<g2f::ParamBaseCRef, std::vector<double>>>(),
                           "loglike"_a, "residuals"_a = std::vector<double>{},
                           "jacobians"_a = std::map<g2f::ParamBaseCRef, std::vector<double>>{})
                      .def_readonly("loglike", &g2f::PriorEvaluation::loglike)
                      .def_readonly("jacobians", &g2f::PriorEvaluation::jacobians)
                      .def_readonly("residuals", &g2f::PriorEvaluation::residuals)
                      .def("__repr__", [](const g2f::PriorEvaluation &self) { return self.repr(true); })
                      .def("__str__", &g2f::PriorEvaluation::str);
    auto _p = py::class_<g2f::Prior, std::shared_ptr<g2f::Prior>, lsst::gauss2d::Object>(m, "Prior");
}
