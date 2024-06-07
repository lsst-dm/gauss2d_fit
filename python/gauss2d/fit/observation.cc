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

#include "lsst/gauss2d/fit/observation.h"
#include "lsst/gauss2d/fit/parametric.h"
#include "lsst/gauss2d/python/pyimage.h"
#include "pybind11.h"

namespace py = pybind11;
using namespace pybind11::literals;

namespace g2f = lsst::gauss2d::fit;
namespace g2p = gauss2d::python;

typedef g2p::PyImage<double> Image;
typedef g2p::PyImage<bool> Mask;
typedef g2f::Observation<double, Image, Mask> Observation;

void bind_observation(py::module &m) {
    auto _o = py::class_<Observation, std::shared_ptr<Observation>, g2f::Parametric>(m, "Observation")
                      .def(py::init<std::shared_ptr<Image>, std::shared_ptr<Image>, std::shared_ptr<Mask>,
                                    const g2f::Channel &>(),
                           "image"_a, "sigma_inv"_a, "mask_inv"_a, "channel"_a)
                      .def_property_readonly("channel", &Observation::get_channel)
                      .def_property_readonly("image", &Observation::get_image)
                      .def_property_readonly("mask_inv", &Observation::get_mask_inverse)
                      .def_property_readonly("sigma_inv", &Observation::get_sigma_inverse)
                      .def_property_readonly("n_cols", &Observation::get_n_cols)
                      .def_property_readonly("n_rows", &Observation::get_n_rows)
                      .def("parameters", &Observation::get_parameters, "parameters"_a = g2f::ParamRefs(),
                           "paramfilter"_a = nullptr)
                      .def("__repr__", [](const Observation &self) { return self.repr(true); })
                      .def("__str__", &Observation::str);
}
