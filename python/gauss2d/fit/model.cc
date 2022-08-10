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

#include "gauss2d/gaussian.h"
#include "gauss2d/fit/model.h"
#include "gauss2d/fit/parametricmodel.h"
#include "gauss2d/python/pyimage.h"

namespace py = pybind11;
using namespace pybind11::literals;

namespace g2f = gauss2d::fit;
namespace g2p = gauss2d::python;

typedef g2f::Model<double, g2p::PyImage<double>, g2p::PyImage<size_t>, g2p::PyImage<bool>> Model;

void bind_model(py::module &m)
{
   auto _o = py::class_<Model, std::shared_ptr<Model>, g2f::ParametricModel
    >(m, "Model")
        .def(py::init<
            std::shared_ptr<const Model::ModelData>,
            Model::PsfModels &,
            Model::Sources &
        >(), "data"_a, "psfmodels"_a, "sources"_a)
        .def("evaluate", &Model::evaluate)
        .def_property_readonly("data", &Model::get_data)
        .def("gaussians", [](const Model & m, const g2f::Channel & c)
            { return std::shared_ptr<const gauss2d::Gaussians>(m.get_gaussians(c)); })
        .def_property_readonly("outputs", &Model::get_outputs)
        .def("parameters", &Model::get_parameters, "parameters"_a=g2f::ParamRefs(), "paramfilter"_a=nullptr)
        .def_property_readonly("psfmodels", &Model::get_psfmodels)
        .def_property_readonly("sources", &Model::get_sources)
        .def("setup_evaluators", &Model::setup_evaluators, "save_gradients"_a=false, "print"_a=false)
        .def("__repr__", &Model::str)
    ;
}
