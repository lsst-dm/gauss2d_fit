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

#include "lsst/gauss2d/fit/sersicmix.h"

#include "pybind11.h"

namespace py = pybind11;
using namespace pybind11::literals;

namespace g2f = lsst::gauss2d::fit;

void bind_sersicmix(py::module &m) {
    auto _is = py::class_<g2f::IntegralSize, std::shared_ptr<g2f::IntegralSize>, lsst::gauss2d::Object>(
                       m, "IntegralSize")
                       .def(py::init<const double, const double>(), "integral"_a = 0, "sigma"_a = 0)
                       .def_readonly("integral", &g2f::IntegralSize::integral)
                       .def_readonly("sigma", &g2f::IntegralSize::sigma)
                       .def("__repr__", [](const g2f::IntegralSize &self) { return self.repr(true); })
                       .def("__str__", &g2f::IntegralSize::str);

    auto _smi = py::class_<g2f::SersicMixInterpolator, std::shared_ptr<g2f::SersicMixInterpolator>>(
            m, "SersicMixInterpolator");

    auto _smv = py::class_<g2f::SersicMixValues, std::shared_ptr<g2f::SersicMixValues>, lsst::gauss2d::Object>(
                        m, "SersicMixValues")
                        .def(py::init<double, std::vector<g2f::IntegralSize>>(), "sersicindex"_a, "values"_a)
                        .def_readonly("sersicindex", &g2f::SersicMixValues::sersicindex)
                        .def_readonly("values", &g2f::SersicMixValues::values)
                        .def("__repr__", [](const g2f::SersicMixValues &self) { return self.repr(true); })
                        .def("__str__", &g2f::SersicMixValues::str);

    m.def("sersic_mix_knots", &g2f::get_sersic_mix_knots_copy, "order"_a = g2f::SERSICMIX_ORDER_DEFAULT,
          py::return_value_policy::copy);
}
