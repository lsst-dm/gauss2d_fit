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

#include "lsst/gauss2d/fit/data.h"
#include "lsst/gauss2d/fit/observation.h"
#include "lsst/gauss2d/fit/parametric.h"
#include "lsst/gauss2d/python/image.h"
#include "pybind11.h"

namespace py = pybind11;
using namespace pybind11::literals;

namespace g2f = lsst::gauss2d::fit;
namespace g2p = lsst::gauss2d::python;

template <typename T>
void declare_data(py::module &m, std::string str_type) {
    typedef g2f::Data<T, g2p::Image<T>, g2p::Image<bool>> Data;
    typedef g2f::Observation<T, g2p::Image<T>, g2p::Image<bool>> Observation;
    std::string pyclass_name = std::string("Data") + str_type;
    py::class_<Data, std::shared_ptr<Data>, g2f::Chromatic, g2f::Parametric>(m, pyclass_name.c_str())
            .def(py::init<std::vector<std::shared_ptr<const Observation>>>(), "data"_a)
            .def_property_readonly("channels", &Data::get_channels)
            .def("parameters", &Data::get_parameters, "parameters"_a = g2f::ParamRefs(),
                 "paramfilter"_a = nullptr)
            .def_property_readonly("size", &Data::size)
            .def("__getitem__", &Data::at)
            .def("__len__", &Data::size)
            //                      .def("__repr__", [](const Data &self) { return self.repr(true);
            //                      })
            .def("__str__", [str_type](const Data &self) {
                std::string str = self.str();
                str = g2p::replace_images_types<T, lsst::gauss2d::idx_type>(str, str_type,
                                                                            self.CC_NAMESPACE_SEPARATOR);
                return str;
            });
}

void bind_data(py::module &m) {
    declare_data<double>(m, "D");
    declare_data<float>(m, "F");
}
