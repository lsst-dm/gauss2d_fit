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

#include "lsst/gauss2d/fit/channel.h"
#include "pybind11.h"

namespace py = pybind11;
using namespace pybind11::literals;

namespace g2f = lsst::gauss2d::fit;

void bind_channel(py::module &m) {
    auto _c = py::class_<g2f::Channel, std::shared_ptr<g2f::Channel> >(m, "Channel")
                      .def(py::init(&g2f::Channel::make))
                      .def_static("erase", &g2f::Channel::erase)
                      .def_property_readonly_static("all",
                                                    [](py::object) { return g2f::Channel::get_channels(); })
                      .def_static("find", &g2f::Channel::find_channel)
                      .def_static("get", &g2f::Channel::get_channel)
                      .def_readonly("name", &g2f::Channel::name)
                      .def_property_readonly_static("NONE",
                                                    [](py::object) { return g2f::Channel::NONE_PTR(); })
                      .def("__repr__", [](const g2f::Channel &self) { return self.repr(true); })
                      .def("__str__", &g2f::Channel::str);
    /*
        const bool operator < ( const Channel &c ) const;
        const bool operator == ( const Channel &c ) const;
        const bool operator != ( const Channel &c ) const;
    */
}
