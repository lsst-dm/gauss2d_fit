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
#include <string>

#include "gauss2d/ellipse.h"
#include "gauss2d/fit/ellipseparameters.h"
#include "gauss2d/fit/parameters.h"
#include "gauss2d/fit/parametric.h"
#include "pybind11.h"

namespace py = pybind11;
using namespace pybind11::literals;

namespace g2f = gauss2d::fit;

void bind_ellipseparameters(py::module &m)
{
    auto _e = py::class_<g2f::EllipseParameters,
        std::shared_ptr<g2f::EllipseParameters>,
        gauss2d::EllipseData, g2f::Parametric
    >(m, "EllipseParameters")
        .def(py::init<double, double, double>(), "sigma_x"_a=0, "sigma_y"_a=0, "rho"_a=0)
        .def(py::init<
            std::shared_ptr<gauss2d::fit::SigmaXParameter>,
            std::shared_ptr<gauss2d::fit::SigmaYParameter>,
            std::shared_ptr<gauss2d::fit::RhoParameter>
            >(),
            "sigma_x"_a=nullptr, "sigma_y"_a=nullptr, "rho"_a=nullptr
        )
        .def_property("rho", &g2f::EllipseParameters::get_rho,
            &g2f::EllipseParameters::set_rho)
        .def_property("sigma_x", &g2f::EllipseParameters::get_sigma_x,
            &g2f::EllipseParameters::set_sigma_x)
        .def_property("sigma_y", &g2f::EllipseParameters::get_sigma_y,
            &g2f::EllipseParameters::set_sigma_y)
        .def_property("xyr", &g2f::EllipseParameters::get_xyr,
            &g2f::EllipseParameters::set_xyr)
        .def_property_readonly("rho_param", &g2f::EllipseParameters::get_rho_param)
        .def_property_readonly("sigma_x_param", &g2f::EllipseParameters::get_sigma_x_param)
        .def_property_readonly("sigma_y_param", &g2f::EllipseParameters::get_sigma_y_param)
        .def_property_readonly("rho_param_ptr", &g2f::EllipseParameters::get_rho_param_ptr)
        .def_property_readonly("sigma_x_param_ptr", &g2f::EllipseParameters::get_sigma_x_param_ptr)
        .def_property_readonly("sigma_y_param_ptr", &g2f::EllipseParameters::get_sigma_y_param_ptr)
        .def("set", &g2f::EllipseParameters::set)
        .def("__repr__", &g2f::EllipseParameters::str)
    ;
}
