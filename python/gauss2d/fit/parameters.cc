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

#include "pybind11.h"
#include "utils.h"

namespace py = pybind11;
using namespace pybind11::literals;

namespace g2f = gauss2d::fit;

template<typename T>
constexpr std::string_view limits_name();

template<> constexpr std::string_view limits_name<double>() { return ", gauss2d.fit.LimitsD";}

// This class exists solely to hold a string so that it can be
// initialized first in an inheritance hierarchy.
// TODO: Find a better solution. Consider private inheritance from std::string.
class Name
{
public:
    const std::string value;
    
    Name(const std::string name) : value(name) {};
};

template<typename T>
class Limits : public Name, public parameters::Limits<T>
{
public:
    Limits(T min = -std::numeric_limits<T>::infinity(), T max = -std::numeric_limits<T>::infinity(),
        const std::string name="") : Name(name), 
            parameters::Limits<T>(min, max, std::string_view(value), limits_name<T>()) {
        };
};

template<typename T>
void declare_limits(py::module &m) {
    using Base = parameters::Limits<T>;
    using Class = Limits<T>;
    std::string pyclass_name = std::string(limits_name<T>().substr(14));
    auto _n = py::class_<Name, std::shared_ptr<Name>>(m, "Name");
    py::class_<Base, std::shared_ptr<Base>>(m, ("_" + pyclass_name).c_str())
    .def("check", &Base::check)
    .def("clip", &Base::clip)
    .def_property("min", &Base::get_min, &Base::set_min)
    .def_property("max", &Base::get_max, &Base::set_max)
    .def("__repr__", &Base::str);

    py::class_<Class, std::shared_ptr<Class>, Name, Base>(m, pyclass_name.c_str())
    .def(py::init<T, T, const std::string>(),
         "min"_a=-std::numeric_limits<T>::infinity(),
         "max"_a=std::numeric_limits<T>::infinity(),
         "name"_a=""
    )
    .def("check", &Class::check)
    .def("clip", &Class::clip)
    .def_property("min", &Class::get_min, &Class::set_min)
    .def_property("max", &Class::get_max, &Class::set_max)
    .def("__repr__", &Class::str);
}

template<typename T, class C>
auto declare_parameter(py::module &m, std::string name, std::string suffix=g2f::suffix_type_str<T>()) {
    using Base = parameters::ParameterBase<T>;
    using Class = parameters::Parameter<T, C>;
    using SetC = typename C::SetC;
    std::string pyclass_name = name + "Parameter" + g2f::suffix_type_str<T>();
    return py::class_<C, std::shared_ptr<C>, Base>(m, pyclass_name.c_str())
    .def(py::init<
            T, std::shared_ptr<const Limits<T>>, std::shared_ptr<const parameters::Transform<T>>,
            std::shared_ptr<const parameters::Unit>, bool, std::string,
            const SetC, const SetC
        >(),
        "value"_a=Class::_get_default(), "limits"_a=nullptr, "transform"_a=nullptr, "unit"_a=gauss2d::fit::unit_none, "fixed"_a=false, "label"_a="",
        // TODO: Can't seem to get nullptr (static_cast or otherwise)/None or equivalent to work here
        "inheritors"_a=SetC(), "modifiers"_a=SetC()
    )
    .def_property_readonly_static("default", [](py::object) { return Class::_get_default(); })
    .def_property_readonly_static("desc", [](py::object) { return Class::_get_desc(); })
    .def_property("fixed", &Class::get_fixed, &Class::set_fixed)
    .def_property("free", &Class::get_free, &Class::set_free)
    .def_property("inheritors", &Class::get_inheritors, &Class::set_inheritors)
    .def_property("label", &Class::get_label, &Class::set_label)
    .def_property("limits", &Class::get_limits, &Class::set_limits, py::keep_alive<1, 2>())
    .def_property_readonly("linear", &Class::get_linear)
    .def_property("modifiers", &Class::get_modifiers, &Class::set_modifiers)
    .def_property_readonly_static("min", [](py::object) { return Class::_get_min(); })
    .def_property_readonly_static("max", [](py::object) { return Class::_get_max(); })
    .def_property_readonly_static("name", [](py::object) { return Class::_get_name(); } )
    .def_property_readonly("ptr", &Class::ptr)
    // TODO: Figure out if it's possible to bind these
    //.def_property_readonly_static("limits_maximal", &Class::limits_maximal)
    .def_property("transform", &Class::get_transform, &Class::set_transform)
    .def_property_readonly("transform_derivative", &Class::get_transform_derivative)
    //.def_property_readonly_static("transform_none", &Class::transform_none)
    .def_property("value", &Class::get_value, &Class::set_value)
    .def_property("value_transformed", &Class::get_value_transformed, &Class::set_value_transformed)

/*
    Additional methods that might be worth wrapping.

    inline std::shared_ptr<Limits<T>> & _get_limits() { return _limits_ptr; }
    inline std::shared_ptr<Transform<T>> & _get_transform() { return _transform_ptr; }
*/
//    This doesn't work - should it?
//    .def(py::self == py::self)
//    .def(py::self != py::self)
//    .def(hash(py::self))
    .def("__repr__", &Class::str);
}

template<typename T>
void declare_transform_base(py::module &m, std::string suffix=g2f::suffix_type_str<T>()) {
    using Class = parameters::Transform<T>;
    py::class_<Class, std::shared_ptr<Class>>(m, ("Transform" + suffix).c_str());
}

template<typename T, class C, bool has_factor, bool has_limits, typename ...Arguments>
void declare_transform_full(
    py::module &m, std::string name, std::string suffix=g2f::suffix_type_str<T>()
) {
    using Class = C;
    auto x = py::class_<Class, std::shared_ptr<Class>, parameters::Transform<T>>(
        m, (name + "Transform" + suffix).c_str()
    )
    .def("description", &Class::description)
    .def("derivative", &Class::derivative)
    .def("forward", &Class::forward)
    .def("reverse", &Class::reverse)
    .def("__repr__", &Class::str);
    if constexpr(has_factor) x.def_property("factor", &Class::get_factor, &Class::set_factor);
    if constexpr(has_limits) x.def_property("limits", &Class::get_limits, &Class::set_limits);
    // TODO: Figure out a neater way to do this
    constexpr const bool has_both = has_factor && has_limits;
    if constexpr(has_both) x.def(py::init<Arguments...>(), "limits"_a=nullptr, "factor"_a=1.);
    else if constexpr(has_factor) x.def(py::init<Arguments...>(), "factor"_a=1.);
    else if constexpr(has_limits) x.def(py::init<Arguments...>(), "limits"_a=nullptr);
    else if constexpr(!has_both) x.def(py::init<>());
}

template<typename T, class C, typename ...Arguments>
void declare_transform(
    py::module &m, std::string name, std::string suffix=g2f::suffix_type_str<T>()
) {
    declare_transform_full<T, C, false, false, Arguments...>(m, name, suffix);
}

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
    declare_parameter<double, gauss2d::fit::RadiusScaleParameter>(m, "RadiusScale");
    declare_parameter<double, gauss2d::fit::RhoParameter>(m, "Rho");
    declare_parameter<double, gauss2d::fit::SersicIndexParameter>(m, "SersicIndex");
    declare_parameter<double, gauss2d::fit::SigmaXParameter>(m, "SigmaX");  
    declare_parameter<double, gauss2d::fit::SigmaYParameter>(m, "SigmaY");
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
