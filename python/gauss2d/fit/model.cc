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

typedef g2p::PyImage<double> Image;
typedef g2f::Model<double, Image, g2p::PyImage<size_t>, g2p::PyImage<bool>> Model;

void bind_model(py::module &m) {
    auto model = py::class_<Model, std::shared_ptr<Model>, g2f::ParametricModel>(m, "Model");
    auto _e = py::enum_<Model::EvaluatorMode>(model, "EvaluatorMode")
                      .value("image", Model::EvaluatorMode::image)
                      .value("loglike", Model::EvaluatorMode::loglike)
                      .value("loglike_image", Model::EvaluatorMode::loglike_image)
                      .value("loglike_grad", Model::EvaluatorMode::loglike_grad)
                      .value("jacobian", Model::EvaluatorMode::jacobian)
                      .export_values();
    model.def(
        py::init<std::shared_ptr<const Model::ModelData>, Model::PsfModels &, Model::Sources &, Model::Priors &>(),
            "data"_a, "psfmodels"_a, "sources"_a, "priors"_a=Model::Priors{})
        .def("evaluate", py::overload_cast<>(&Model::evaluate))
        .def_property_readonly("data", &Model::get_data)
        .def("gaussians",
             [](const Model &m, const g2f::Channel &c) {
                 return std::shared_ptr<const gauss2d::Gaussians>(m.get_gaussians(c));
             })
        .def_property_readonly("mode", &Model::get_mode)
        .def_property_readonly("outputs", &Model::get_outputs)
        .def("parameters", &Model::get_parameters, "parameters"_a = g2f::ParamRefs(),
             "paramfilter"_a = nullptr)
        .def_property_readonly("priors", &Model::get_priors)
        .def_property_readonly("psfmodels", &Model::get_psfmodels)
        .def_property_readonly("sources", &Model::get_sources)
        .def("setup_evaluators", &Model::setup_evaluators,
             "evaluatormode"_a = Model::EvaluatorMode::image,
             "outputs"_a = std::vector<std::vector<std::shared_ptr<Image>>>{},
             "residuals"_a = std::vector<std::shared_ptr<Image>>{},
             "outputs_prior"_a = std::vector<std::shared_ptr<Image>>{},
             "residuals_prior"_a = std::shared_ptr<Image>{},
             "force"_a = false,
             "print"_a = false)
        .def("verify_jacobian", &Model::verify_jacobian, "findiff_frac"_a = 1e-4, "findiff_add"_a = 1e-4,
             "rtol"_a = 1e-3, "atol"_a = 1e-3)
        .def("__repr__", [](const Model &self) { return self.repr(true); })
        .def("__str__", &Model::str);
}
