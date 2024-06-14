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

#ifndef LSST_GAUSS2D_FIT_PYTHON_PARAMETERS_H
#define LSST_GAUSS2D_FIT_PYTHON_PARAMETERS_H

#include <pybind11/attr.h>
#include <pybind11/pybind11.h>
#include <pybind11/operators.h>
#include <pybind11/stl.h>

#include <memory>
#include <string>

#include "lsst/modelfit/parameters.h"

#include "gauss2d/fit/parameters.h"
#include "gauss2d/fit/transforms.h"

#include "pybind11.h"
#include "utils.h"

namespace py = pybind11;
using namespace pybind11::literals;

namespace g2f = lsst::gauss2d::fit;
namespace parameters = lsst::modelfit::parameters;

template <typename T>
void declare_limits(py::module &m) {
    using Class = parameters::Limits<T>;
    std::string pyclass_name = std::string("Limits") + g2f::suffix_type_str<T>();
    py::class_<Class, std::shared_ptr<Class>>(m, pyclass_name.c_str())
            .def(py::init<T, T, const std::string>(), "min"_a = -std::numeric_limits<T>::infinity(),
                 "max"_a = std::numeric_limits<T>::infinity(), "name"_a = "")
            .def("check", &Class::check)
            .def("clip", &Class::clip)
            .def_property("min", &Class::get_min, &Class::set_min)
            .def_property("max", &Class::get_max, &Class::set_max)
            .def("__repr__", [](const Class &self) { return self.repr(true); })
            .def("__str__", &Class::str);
}

/*
    Class is the CRTP class for CRTP parameters, whereas C is the "derived" class.
    For classes derived from a CRTP parameter, the two are the same.
    TODO: Come up with a better naming scheme for this.
*/
template <class Class, class C, typename... Args>
auto declare_parameter_methods(py::class_<C, Args...> c) {
    return c.def_property_readonly("default", &Class::get_default)
            .def_property_readonly("desc", &Class::get_desc)
            .def_property("fixed", &Class::get_fixed, &Class::set_fixed)
            .def_property("free", &Class::get_free, &Class::set_free)
            .def_property("label", &Class::get_label, &Class::set_label)
            // Return a copy of the limits, because the C++ func returns a const Limits &
            // and calling setters on it in Python could cause segfaults
            .def_property(
                    "limits",
                    [](const C &self) {
                        return parameters::Limits<double>{self.get_limits().get_min(),
                                                          self.get_limits().get_max()};
                    },
                    &Class::set_limits)
            .def_property_readonly("linear", &Class::get_linear)
            .def_property_readonly("min", &Class::get_min)
            .def_property_readonly("max", &Class::get_max)
            .def_property_readonly("name", &Class::get_name)
            .def_property_readonly("ptr", &Class::ptr)
            // TODO: Figure out if it's possible to bind these
            //.def_property_readonly_static("limits_maximal", &Class::limits_maximal)
            // TODO: Determine if this will also segfault when mutating the transform, since
            // the C++ returns a const ref (most but not all transforms are immutable)
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
            .def("__str__", &Class::str);
}

template <typename T, class C, class... Bases>
auto declare_parameter(py::module &m, std::string name, std::string suffix = g2f::suffix_type_str<T>()) {
    using Base = parameters::ParameterBase<T>;
    using Class = parameters::Parameter<T, C>;
    return declare_parameter_methods<Class, C, std::shared_ptr<C>, Base>(
            py::class_<C, std::shared_ptr<C>, Base, Bases...>(m, (name + "Parameter" + suffix).c_str())
                    .def(py::init<T, std::shared_ptr<const parameters::Limits<T>>,
                                  std::shared_ptr<const parameters::Transform<T>>,
                                  std::shared_ptr<const parameters::Unit>, bool, std::string>(),
                         "value"_a = Class::_get_default(), "limits"_a = nullptr, "transform"_a = nullptr,
                         "unit"_a = lsst::gauss2d::fit::unit_none, "fixed"_a = false, "label"_a = "")
                    .def("__repr__", [](const C &self) { return self.repr(true, "."); }));
}

template <typename T, class C, class... Bases>
auto declare_sizeparameter(py::module &m, std::string name) {
    return declare_parameter<T, C, Bases...>(m, name).def_property("size", &C::get_size, &C::set_size);
}

template <typename T, class ClassX, class ClassY>
auto declare_sizeparameter_base(py::module &m, std::string suffix = g2f::suffix_type_str<T>()) {
    py::class_<ClassX, std::shared_ptr<ClassX>>(m, ("SizeXParameter" + suffix).c_str());
    py::class_<ClassY, std::shared_ptr<ClassY>>(m, ("SizeYParameter" + suffix).c_str());
}

template <typename T>
void declare_transform_base(py::module &m) {
    using Class = parameters::Transform<T>;
    py::class_<Class, std::shared_ptr<Class>>(m, "TransformD");
}

template <typename T, class C, bool has_factor, bool has_limits, typename... Arguments>
void declare_transform_full(py::module &m, std::string name) {
    using Class = C;
    auto x = py::class_<Class, std::shared_ptr<Class>, parameters::Transform<T>>(
                     m, (name + "TransformD").c_str())
                     .def("description", &Class::description)
                     .def("derivative", &Class::derivative)
                     .def("forward", &Class::forward)
                     .def("reverse", &Class::reverse)
                     .def("__repr__", [](const Class &self) { return self.repr(true, "."); })
                     .def("__str__", &Class::str);
    if constexpr (has_factor) x.def_property("factor", &Class::get_factor, &Class::set_factor);
    if constexpr (has_limits) x.def_property("limits", &Class::get_limits, &Class::set_limits);
    // TODO: Figure out a neater way to do this
    constexpr const bool has_both = has_factor && has_limits;
    if constexpr (has_both)
        x.def(py::init<Arguments...>(), "limits"_a = nullptr, "factor"_a = 1.);
    else if constexpr (has_factor)
        x.def(py::init<Arguments...>(), "factor"_a = 1.);
    else if constexpr (has_limits)
        x.def(py::init<Arguments...>(), "limits"_a = nullptr);
    else if constexpr (!has_both)
        x.def(py::init<>());
}

template <typename T, class C, typename... Arguments>
void declare_transform(py::module &m, std::string name) {
    declare_transform_full<T, C, false, false, Arguments...>(m, name);
}

#endif
