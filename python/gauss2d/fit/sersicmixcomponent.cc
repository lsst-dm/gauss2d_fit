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

#include "gauss2d/fit/component.h"
#include "gauss2d/fit/parameters.h"
#include "gauss2d/fit/sersicmixcomponent.h"
#include "gauss2d/fit/sersicparametricellipse.h"

#include "parameters.h"
#include "parameters/parameter.h"
#include "pybind11.h"

namespace py = pybind11;
using namespace pybind11::literals;

namespace g2f = gauss2d::fit;

void bind_sersicmixcomponent(py::module &m)
{
    using T = double;
    using C = g2f::SersicMixComponentIndexParameter;
    using Base = g2f::SersicIndexParameter;
    using SetC = typename g2f::SersicMixComponentIndexParameter::SetC;

    std::string pyclass_name = "SersicMixComponentIndexParameter" + g2f::suffix_type_str<T>();
    declare_parameter_methods<C, C, std::shared_ptr<C>, Base>(
        // note that Base is the actual name of the base Parameter class, not the CRTP class
        // that it is "derived" from
        py::class_<C, std::shared_ptr<C>, Base>(m, pyclass_name.c_str())
    )
    // new properties
        .def_property_readonly("integralratio", &g2f::SersicMixComponentIndexParameter::get_integralratio)
        .def_readonly("order", &g2f::SersicMixComponentIndexParameter::order)
        .def_property_readonly("sizeratio", &g2f::SersicMixComponentIndexParameter::get_sizeratio)
    // constructor with added arg
        .def(py::init<
            T, std::shared_ptr<const Limits<T>>, std::shared_ptr<const parameters::Transform<T>>,
            std::shared_ptr<const parameters::Unit>, bool, std::string,
            const SetC, const SetC, const std::shared_ptr<const g2f::SersicMixInterpolator>
        >(),
        "value"_a=g2f::SersicMixComponentIndexParameter::_get_default(),
        "limits"_a=nullptr, "transform"_a=nullptr, "unit"_a=gauss2d::fit::unit_none,
        "fixed"_a=false, "label"_a="",
        // TODO: Can't seem to get nullptr (static_cast or otherwise)/None or equivalent to work here
        "inheritors"_a=SetC(), "modifiers"_a=SetC(),
        "interpolator"_a=nullptr
    );

    auto _e = py::class_<g2f::SersicMixComponent,
        std::shared_ptr<g2f::SersicMixComponent>,
        g2f::Component
    >(m, "SersicMixComponent")
        .def(py::init<
                std::shared_ptr<g2f::SersicParametricEllipse>,
                std::shared_ptr<g2f::CentroidParameters>,
                std::shared_ptr<g2f::IntegralModel>,
                std::shared_ptr<g2f::SersicMixComponentIndexParameter>
            >(),
            "ellipse"_a=nullptr, "centroid"_a=nullptr, "integral"_a=nullptr, "sersicindex"_a=nullptr)
        .def("parameters", &g2f::SersicMixComponent::get_parameters,
            "parameters"_a=g2f::ParamRefs(), "paramfilter"_a=nullptr)
        .def("gaussians", [](const g2f::SersicMixComponent &g, const g2f::Channel & c)
            { return std::shared_ptr<const gauss2d::Gaussians>(g.get_gaussians(c)); })
        .def("__repr__", &g2f::SersicMixComponent::str)
    ;
}
