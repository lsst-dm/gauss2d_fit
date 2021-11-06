/*
 * This file is part of gauss2d.
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

#include <memory>
#include <pybind11/attr.h>
#ifndef GAUSS2D_CENTROID_H
#include "gauss2d/centroid.h"
#endif

#ifndef GAUSS2D_ELLIPSE_H
#include "gauss2d/ellipse.h"
#endif

#ifndef PARAMETERS_TRANSFORM_H
#include "parameters/transform.h"
#endif

#ifndef GAUSS2D_FIT_ELLIPSEPARAMETERS_H
#include "gauss2d/fit/ellipseparameters.h"
#endif

#ifndef GAUSS2D_FIT_PARAMETERS_H
#include "gauss2d/fit/parameters.h"
#endif

#ifndef GAUSS2DFIT_TRANSFORMS_H
#include "gauss2d/fit/transforms.h"
#endif

#include <limits>
//#include <iostream>
#include <string>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;
using namespace pybind11::literals;

// This is so annoying that I'm considering borrowing a solution for static string 
template<typename T>
constexpr std::string_view suffix_type();

template<> constexpr std::string_view suffix_type<double>() { return "D";}

template<typename T>
std::string suffix_type_str() { return std::string(suffix_type<T>()); }

template<typename T>
constexpr std::string_view limits_name();

template<> constexpr std::string_view limits_name<double>() { return ", gauss2d.fit.LimitsD";}

// Yes, this class exists solely to hold a string so that it can be initialized first
// in an inheritance hierarchy. Solve this better if you can.
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
    py::class_<Base, std::shared_ptr<Base>>(m, ("_" + pyclass_name).c_str())
    .def("check", &Base::check)
    .def("clip", &Base::clip)
    .def_property("min", &Base::get_min, &Base::set_min)
    .def_property("max", &Base::get_max, &Base::set_max)
    .def("__repr__", &Base::str);
    py::class_<Class, std::shared_ptr<Class>, Base>(m, pyclass_name.c_str())
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
void declare_parameter(py::module &m, std::string name, std::string suffix=suffix_type_str<T>()) {
    using Class = parameters::Parameter<T, C>;
    std::string pyclass_name = name + "Parameter" + suffix_type_str<T>();
    py::class_<C, std::shared_ptr<C>>(m, pyclass_name.c_str())
    .def(py::init<T, std::shared_ptr<const Limits<T>>, std::shared_ptr<const parameters::Transform<T>>>(),
        "value"_a=Class::get_default(), "limits"_a=nullptr, "transform"_a=nullptr)
    .def_property("limits", &Class::get_limits, &Class::set_limits)
    .def_property_readonly("min", &Class::get_min)
    .def_property_readonly("max", &Class::get_max)
    .def_property_readonly("ptr", &Class::ptr)
    // TODO: Figure out if it's possible to bind these
    //.def_property_readonly("limits_maximal", &Class::limits_maximal)
    .def_property("transform", &Class::get_transform, &Class::set_transform)
    //.def_property_readonly("transform_none", &Class::transform_none)
    .def_property("value", &Class::get_value, &Class::set_value)
    .def_property("value_transformed", &Class::get_value_transformed, &Class::set_value_transformed)

/*
    Additional methods that might be worth wrapping.

    inline std::shared_ptr<Limits<T>> & _get_limits() { return _limits_ptr; }
    inline std::shared_ptr<Transform<T>> & _get_transform() { return _transform_ptr; }
    static const std::string get_desc() { return _desc_t<C>(0); }
    static constexpr T get_min() { return _min_t<C>(0); }
    static constexpr T get_max() { return _max_t<C>(0); }
    static inline const std::string get_name() { return _name_t<C>(0); }
    static inline const std::string get_type_name() { return std::string(type_name<C>()); }

    static constexpr const Limits<T> limits_maximal = Limits<T>(get_min(), get_max(), type_name<C>(), ".limits_maximal");
    static constexpr const UnitTransform<T> transform_none = UnitTransform<T>();
*/
    .def("__repr__", &Class::str);
}

template<typename T>
void declare_transform_base(py::module &m, std::string suffix=suffix_type_str<T>()) {
    using Class = parameters::Transform<T>;
    py::class_<Class, std::shared_ptr<Class>>(m, ("Transform" + suffix).c_str());
}

template<typename T, class C, bool has_factor, bool has_limits, typename ...Arguments>
void declare_transform_full(
    py::module &m, std::string name, std::string suffix=suffix_type_str<T>()
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
    py::module &m, std::string name, std::string suffix=suffix_type_str<T>()
) {
    declare_transform_full<T, C, false, false, Arguments...>(m, name, suffix);
}

PYBIND11_MODULE(_gauss2dfit, m)
{
    m.doc() = "Gauss2D::Fit Pybind11 functions"; // optional module docstring

    py::class_<gauss2d::fit::CentroidParameters,
        std::shared_ptr<gauss2d::fit::CentroidParameters>,
        gauss2d::CentroidData
    >(m, "CentroidParameters")
        .def(py::init<double, double>(), "x"_a=0, "y"_a=0)
        .def(py::init<
            std::shared_ptr<gauss2d::fit::CentroidXParameter>,
            std::shared_ptr<gauss2d::fit::CentroidYParameter>>(),
            "x"_a=nullptr, "y"_a=nullptr
        )
        .def_property("x", &gauss2d::fit::CentroidParameters::get_x, &gauss2d::fit::CentroidParameters::set_x)
        .def_property("y", &gauss2d::fit::CentroidParameters::get_y, &gauss2d::fit::CentroidParameters::set_y)
        .def_property("xy", &gauss2d::fit::CentroidParameters::get_xy, 
            &gauss2d::fit::CentroidParameters::set_xy)
        .def_property_readonly("get_x_param", &gauss2d::fit::CentroidParameters::get_x_param)
        .def_property_readonly("get_y_param", &gauss2d::fit::CentroidParameters::get_y_param)
        .def_property_readonly("get_x_param_ptr", &gauss2d::fit::CentroidParameters::get_x_param_ptr)
        .def_property_readonly("get_y_param_ptr", &gauss2d::fit::CentroidParameters::get_y_param_ptr)
        .def("__repr__", &gauss2d::fit::CentroidParameters::str)
    ;
    py::class_<gauss2d::fit::EllipseParameters,
        std::shared_ptr<gauss2d::fit::EllipseParameters>,
        gauss2d::EllipseData
    >(m, "Ellipse")
        .def(py::init<double, double, double>(), "sigma_x"_a=0, "sigma_y"_a=0, "rho"_a=0)
        .def_property("rho", &gauss2d::fit::EllipseParameters::get_rho, 
            &gauss2d::fit::EllipseParameters::set_rho)
        .def_property("sigma_x", &gauss2d::fit::EllipseParameters::get_sigma_x,
            &gauss2d::fit::EllipseParameters::set_sigma_x)
        .def_property("sigma_y", &gauss2d::fit::EllipseParameters::get_sigma_y,
            &gauss2d::fit::EllipseParameters::set_sigma_y)
        .def_property("xyr", &gauss2d::fit::EllipseParameters::get_xyr,
            &gauss2d::fit::EllipseParameters::set_xyr)
        .def_property_readonly("rho_param", &gauss2d::fit::EllipseParameters::get_rho_param)
        .def_property_readonly("sigma_x_param", &gauss2d::fit::EllipseParameters::get_sigma_x_param)
        .def_property_readonly("sigma_y_param", &gauss2d::fit::EllipseParameters::get_sigma_y_param)
        .def_property_readonly("rho_param_ptr", &gauss2d::fit::EllipseParameters::get_rho_param_ptr)
        .def_property_readonly("sigma_x_param_ptr", &gauss2d::fit::EllipseParameters::get_sigma_x_param_ptr)
        .def_property_readonly("sigma_y_param_ptr", &gauss2d::fit::EllipseParameters::get_sigma_y_param_ptr)
        .def("set", &gauss2d::fit::EllipseParameters::set)
        .def("__repr__", &gauss2d::fit::EllipseParameters::str)
    ;
    declare_limits<double>(m);
    declare_transform_base<double>(m);
    declare_transform<double, parameters::UnitTransform<double>>(m, "Unit");
    declare_transform<double, gauss2d::fit::InverseTransform>(m, "Inverse");
    declare_transform<double, gauss2d::fit::LogTransform>(m, "Log");
    declare_transform<double, gauss2d::fit::Log10Transform>(m, "Log10");
    declare_transform<double, gauss2d::fit::LogitTransform>(m, "Logit");
    // TODO: Determine why this won't work with std::shared_ptr<parameters::Limits<double>>
    declare_transform_full<double, gauss2d::fit::LogitLimitedTransform, true, true,
        std::shared_ptr<Limits<double>>, double >(m, "LogitLimited");
    declare_parameter<double, gauss2d::fit::IntegralParameter>(m, "Integral");
    declare_parameter<double, gauss2d::fit::CentroidXParameter>(m, "CentroidX");
    declare_parameter<double, gauss2d::fit::CentroidYParameter>(m, "CentroidY");
    declare_parameter<double, gauss2d::fit::RhoParameter>(m, "Rho");
    declare_parameter<double, gauss2d::fit::SigmaXParameter>(m, "SigmaX");
    declare_parameter<double, gauss2d::fit::SigmaYParameter>(m, "SigmaY");
}
