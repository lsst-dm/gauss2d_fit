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
#include <pybind11/operators.h>
#include <pybind11/stl.h>

#include <memory>
#include <string>

#include "gauss2d/fit/parameters.h"
#include "gauss2d/fit/transforms.h"
#include "parameters/transform.h"

#include "parameters.h"
#include "pybind11.h"
#include "utils.h"

namespace py = pybind11;
using namespace pybind11::literals;


void bind_parameters(py::module &m)
{
   auto _u = py::class_<parameters::Unit,
        std::shared_ptr<parameters::Unit>
    >(m, "Unit");
    py::class_<gauss2d::fit::UnitNone,
        std::shared_ptr<gauss2d::fit::UnitNone>,
        parameters::Unit
    >(m, "UnitNone")
        .def(py::init<>())
        .def_property_readonly("name", &gauss2d::fit::UnitNone::get_name)
    ;
    declare_limits<double>(m);
    using Parameter = parameters::ParameterBase<double>;
    auto _p = py::class_<Parameter, std::shared_ptr<Parameter>>(m, "Parameter");
    auto integral = declare_parameter<double, gauss2d::fit::IntegralParameter>(m, "Integral");
    integral.def_property("label", &gauss2d::fit::IntegralParameter::get_label,
        &gauss2d::fit::IntegralParameter::set_label);
    declare_parameter<double, gauss2d::fit::CentroidXParameter>(m, "CentroidX");
    declare_parameter<double, gauss2d::fit::CentroidYParameter>(m, "CentroidY");
    declare_parameter<double, gauss2d::fit::MoffatConcentrationParameter>(m, "MoffatConcentration");
    auto propfrac = declare_parameter<double, gauss2d::fit::ProperFractionParameter>(m, "ProperFraction");
    propfrac.def_property("label", &gauss2d::fit::ProperFractionParameter::get_label,
        &gauss2d::fit::ProperFractionParameter::set_label);
    declare_sizeparameter_base<double>(m);
    declare_parameter<double, gauss2d::fit::RadiusScaleParameter>(m, "RadiusScale");
    declare_sizeparameter<double, gauss2d::fit::ReffXParameter, g2f::SizeXParameter>(m, "ReffX");
    declare_sizeparameter<double, gauss2d::fit::ReffYParameter, g2f::SizeYParameter>(m, "ReffY");
    declare_parameter<double, gauss2d::fit::RhoParameter>(m, "Rho");
    declare_parameter<double, gauss2d::fit::SersicIndexParameter>(m, "SersicIndex");
    declare_sizeparameter<double, gauss2d::fit::SigmaXParameter, g2f::SizeXParameter>(m, "SigmaX"); 
    declare_sizeparameter<double, gauss2d::fit::SigmaYParameter, g2f::SizeYParameter>(m, "SigmaY");
    declare_transform_base<double>(m);
    declare_transform<double, parameters::UnitTransform<double>>(m, "Unit");
    declare_transform<double, gauss2d::fit::InverseTransform>(m, "Inverse");
    declare_transform<double, gauss2d::fit::LogTransform>(m, "Log");
    declare_transform<double, gauss2d::fit::Log10Transform>(m, "Log10");
    declare_transform<double, gauss2d::fit::LogitTransform>(m, "Logit");
    // TODO: Determine why this won't work with std::shared_ptr<parameters::Limits<double>>
    declare_transform_full<double, gauss2d::fit::LogitLimitedTransform, true, true,
        std::shared_ptr<Limits<double>>, double >(m, "LogitLimited");
}
