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

#include "pybind11.h"

#include "lsst/gauss2d/fit/linearsersicmixinterpolator.h"
#include "lsst/gauss2d/fit/sersicmix.h"

namespace py = pybind11;
using namespace pybind11::literals;

namespace g2f = lsst::gauss2d::fit;

void bind_linearsersicmixinterpolator(py::module &m) {
    auto _e = py::class_<g2f::LinearSersicMixInterpolator, std::shared_ptr<g2f::LinearSersicMixInterpolator>,
                         g2f::SersicMixInterpolator>(m, "LinearSersicMixInterpolator")
                      .def(py::init<short>(), "order"_a = g2f::SERSICMIX_ORDER_DEFAULT)
                      .def("integralsizes", &g2f::LinearSersicMixInterpolator::get_integralsizes)
                      .def("integralsizes_derivs",
                           &g2f::LinearSersicMixInterpolator::get_integralsizes_derivs)
                      .def_property_readonly("order", &g2f::LinearSersicMixInterpolator::get_order)
                      .def("__repr__",
                           [](const g2f::LinearSersicMixInterpolator &self) { return self.repr(true); })
                      .def("__str__", &g2f::LinearSersicMixInterpolator::str);
}
