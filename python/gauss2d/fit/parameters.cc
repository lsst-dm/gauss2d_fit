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
#include <vector>

#include "lsst/gauss2d/fit/parameters.h"
#include "lsst/gauss2d/fit/transforms.h"
#include "lsst/gauss2d/fit/util.h"
#include "lsst/modelfit/parameters.h"

#include "parameters.h"
#include "pybind11.h"
#include "utils.h"

namespace py = pybind11;
using namespace pybind11::literals;

namespace g2f = lsst::gauss2d::fit;

void bind_parameters(py::module &m) {
    auto _u = py::class_<parameters::Unit, std::shared_ptr<parameters::Unit>>(m, "Unit");
    py::class_<g2f::UnitNone, std::shared_ptr<g2f::UnitNone>, parameters::Unit>(m, "UnitNone")
            .def(py::init<>())
            .def_property_readonly("name", &g2f::UnitNone::get_name);
    declare_limits<double>(m);
    using Parameter = parameters::ParameterBase<double>;
    auto _p = py::class_<Parameter, std::shared_ptr<Parameter>>(m, "ParameterD");
    auto integral = declare_parameter<double, g2f::IntegralParameterD>(m, "Integral");
    integral.def_property("label", &g2f::IntegralParameterD::get_label, &g2f::IntegralParameterD::set_label);
    declare_parameter<double, g2f::CentroidXParameterD>(m, "CentroidX");
    declare_parameter<double, g2f::CentroidYParameterD>(m, "CentroidY");
    declare_parameter<double, g2f::MeanParameterD>(m, "Mean");
    declare_parameter<double, g2f::MoffatConcentrationParameterD>(m, "MoffatConcentration");
    auto propfrac = declare_parameter<double, g2f::ProperFractionParameterD>(m, "ProperFraction");
    propfrac.def_property("label", &g2f::ProperFractionParameterD::get_label,
                          &g2f::ProperFractionParameterD::set_label);
    declare_sizeparameter_base<double, g2f::SizeXParameterD, g2f::SizeYParameterD>(m);
    declare_parameter<double, g2f::RadiusScaleParameterD>(m, "RadiusScale");
    declare_sizeparameter<double, g2f::ReffXParameterD, g2f::SizeXParameterD>(m, "ReffX");
    declare_sizeparameter<double, g2f::ReffYParameterD, g2f::SizeYParameterD>(m, "ReffY");
    declare_parameter<double, g2f::RhoParameterD>(m, "Rho");
    declare_parameter<double, g2f::SersicIndexParameterD>(m, "SersicIndex");
    declare_sizeparameter<double, g2f::SigmaXParameterD, g2f::SizeXParameterD>(m, "SigmaX");
    declare_sizeparameter<double, g2f::SigmaYParameterD, g2f::SizeYParameterD>(m, "SigmaY");
    declare_parameter<double, g2f::StdDevParameterD>(m, "StdDev");
    declare_transform_base<double>(m);
    declare_transform<double, lsst::modelfit::parameters::UnitTransform<double>>(m, "Unit");
    declare_transform<double, g2f::InverseTransform>(m, "Inverse");
    declare_transform<double, g2f::JanskyToABMagTransform>(m, "JanskyToABMag");
    declare_transform<double, g2f::NanojanskyToABMagTransform>(m, "NanojanskyToABMag");
    declare_transform<double, g2f::LogTransform>(m, "Log");
    declare_transform<double, g2f::Log10Transform>(m, "Log10");
    declare_transform<double, g2f::LogitTransform>(m, "Logit");
    // TODO: Determine why this won't work with std::shared_ptr<parameters::Limits<double>>
    declare_transform_full<double, g2f::LogitLimitedTransform, true, true,
                           std::shared_ptr<parameters::Limits<double>>, double>(m, "LogitLimited");
    const std::vector<std::reference_wrapper<Parameter>> _default = {};
    m.def("params_unique", &g2f::nonconsecutive_unique<std::reference_wrapper<Parameter>>,
          "params"_a = _default);
}
