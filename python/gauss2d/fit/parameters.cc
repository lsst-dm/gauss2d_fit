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

namespace g2f = gauss2d::fit;

void bind_parameters(py::module &m) {
    auto _u = py::class_<parameters::Unit, std::shared_ptr<parameters::Unit>>(m, "Unit");
    py::class_<g2f::UnitNone, std::shared_ptr<g2f::UnitNone>, parameters::Unit>(m, "UnitNone")
            .def(py::init<>())
            .def_property_readonly("name", &g2f::UnitNone::get_name);
    declare_limits<double>(m);
    using Parameter = parameters::ParameterBase<double>;
    auto _p = py::class_<Parameter, std::shared_ptr<Parameter>>(m, "Parameter");
    auto integral = declare_parameter<double, g2f::IntegralParameter>(m, "Integral");
    integral.def_property("label", &g2f::IntegralParameter::get_label, &g2f::IntegralParameter::set_label);
    declare_parameter<double, g2f::CentroidXParameter>(m, "CentroidX");
    declare_parameter<double, g2f::CentroidYParameter>(m, "CentroidY");
    declare_parameter<double, g2f::MeanParameter>(m, "Mean");
    declare_parameter<double, g2f::MoffatConcentrationParameter>(m, "MoffatConcentration");
    auto propfrac = declare_parameter<double, g2f::ProperFractionParameter>(m, "ProperFraction");
    propfrac.def_property("label", &g2f::ProperFractionParameter::get_label,
                          &g2f::ProperFractionParameter::set_label);
    declare_sizeparameter_base<double>(m);
    declare_parameter<double, g2f::RadiusScaleParameter>(m, "RadiusScale");
    declare_sizeparameter<double, g2f::ReffXParameter, g2f::SizeXParameter>(m, "ReffX");
    declare_sizeparameter<double, g2f::ReffYParameter, g2f::SizeYParameter>(m, "ReffY");
    declare_parameter<double, g2f::RhoParameter>(m, "Rho");
    declare_parameter<double, g2f::SersicIndexParameter>(m, "SersicIndex");
    declare_sizeparameter<double, g2f::SigmaXParameter, g2f::SizeXParameter>(m, "SigmaX");
    declare_sizeparameter<double, g2f::SigmaYParameter, g2f::SizeYParameter>(m, "SigmaY");
    declare_parameter<double, g2f::StdDevParameter>(m, "StdDev");
    declare_transform_base<double>(m);
    declare_transform<double, parameters::UnitTransform<double>>(m, "Unit");
    declare_transform<double, g2f::InverseTransform>(m, "Inverse");
    declare_transform<double, g2f::JanskyToABMagTransform>(m, "JanskyToABMag");
    declare_transform<double, g2f::NanojanskyToABMagTransform>(m, "NanojanskyToABMag");
    declare_transform<double, g2f::LogTransform>(m, "Log");
    declare_transform<double, g2f::Log10Transform>(m, "Log10");
    declare_transform<double, g2f::LogitTransform>(m, "Logit");
    // TODO: Determine why this won't work with std::shared_ptr<parameters::Limits<double>>
    declare_transform_full<double, g2f::LogitLimitedTransform, true, true,
                           std::shared_ptr<parameters::Limits<double>>, double>(m, "LogitLimited");
}
