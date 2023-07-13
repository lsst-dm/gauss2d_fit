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

#include "gauss2d/object.h"
#include "gauss2d/fit/param_defs.h"
#include "gauss2d/fit/prior.h"
#include "pybind11.h"

namespace py = pybind11;
using namespace pybind11::literals;

namespace g2f = gauss2d::fit;

void bind_prior(py::module &m) {
    auto _e = py::class_<g2f::PriorEvaluation, std::shared_ptr<g2f::PriorEvaluation>>(m, "PriorEvaluation")
                      .def(py::init<double, std::map<g2f::ParamBaseCRef, std::vector<double>>,
                                    std::vector<double>>(),
                           "loglike"_a, "jacobians"_a = std::map<g2f::ParamBaseCRef, std::vector<double>>{},
                           "residuals"_a = std::vector<double>{})
                      .def_readwrite("loglike", &g2f::PriorEvaluation::loglike)
                      .def_readwrite("jacobians", &g2f::PriorEvaluation::jacobians)
                      .def_readwrite("residuals", &g2f::PriorEvaluation::residuals);
    auto _p = py::class_<g2f::Prior, std::shared_ptr<g2f::Prior>, gauss2d::Object>(m, "Prior");
}
