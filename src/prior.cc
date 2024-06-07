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

#include <string>

#include "lsst/gauss2d/to_string.h"
#include "lsst/gauss2d/type_name.h"

#include "lsst/gauss2d/fit/param_defs.h"
#include "lsst/gauss2d/fit/prior.h"
#include "lsst/gauss2d/fit/util.h"

namespace lsst::gauss2d::fit {

std::string str_jacobians(const std::map<ParamBaseCRef, std::vector<double>>& jacobians,
                          bool python_style = true) {
    std::string str = "{";
    const std::string prefix = python_style ? "" : "{ ";
    const std::string separator = python_style ? ": " : ", ";
    for (const auto& [obj, value] : jacobians) {
        str += obj.get().str() + separator + to_string_float_iter(value) + ", ";
    }
    return str.substr(0, str.size() - 2 * (jacobians.size() > 0)) + "}";
}

PriorEvaluation::PriorEvaluation(double loglike_, std::vector<double> residuals_,
                                 PriorEvaluation::Jacobians jacobians_, bool check_size)
        : loglike(loglike_), residuals(residuals_), jacobians(jacobians_) {
    const auto n_resid = residuals.size();
    if (!((n_resid > 0) || (jacobians.size() == 0))) {
        throw std::invalid_argument("residuals cannot be empty unless jacobians are too");
    }
    if (check_size) {
        for (const auto& [paramref, values] : jacobians) {
            if (values.size() != n_resid) {
                throw std::invalid_argument(paramref.get().str()
                                            + " jacobian size=" + std::to_string(values.size())
                                            + " != n_resid=" + std::to_string(n_resid));
            }
        }
    } else {
    }
}

PriorEvaluation::~PriorEvaluation(){};

double PriorEvaluation::compute_dloglike_dx(const ParamBase& param, bool transformed) const {
    double dll_dx = 0;
    auto jacobian = this->jacobians.find(param);
    if (jacobian != this->jacobians.end()) {
        // jacobian is dresidual/dx
        // loglike is -residual^2/2
        // dloglike/dx = -residual*dresidual/dx
        const auto& jacobian_values = (*jacobian).second;
        for (size_t idx = 0; idx < jacobian_values.size(); ++idx) {
            dll_dx -= jacobian_values[idx] * residuals.at(idx);
        }
    }
    if (!transformed) dll_dx /= param.get_transform_derivative();
    return dll_dx;
}

std::string PriorEvaluation::repr(bool name_keywords, std::string_view namespace_separator) const {
    return type_name_str<PriorEvaluation>(false, namespace_separator) + "("
           + (name_keywords ? "loglike=" : "") + to_string_float(loglike) + ", "
           + (name_keywords ? "residuals=" : "") + to_string_float_iter(residuals) + ", "
           + (name_keywords ? "jacobians=" : "") + str_jacobians(jacobians) + ")";
}

std::string PriorEvaluation::str() const {
    return type_name_str<PriorEvaluation>(true) + "(loglike=" + to_string_float(loglike) + ", residuals="
           + to_string_float_iter(residuals) + ", jacobians=" + str_jacobians(jacobians) + ")";
}

}  // namespace lsst::gauss2d::fit
