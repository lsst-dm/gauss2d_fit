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

#include <cmath>
#include <string>

#include "lsst/gauss2d/to_string.h"
#include "lsst/gauss2d/type_name.h"

#include "lsst/gauss2d/fit/gaussianprior.h"

namespace lsst::gauss2d::fit {
GaussianPrior::GaussianPrior(std::shared_ptr<const ParamBase> param, double mean, double stddev,
                             bool transformed)
        : _param(std::move(param)), _transformed(transformed) {
    if (_param == nullptr) throw std::invalid_argument("param must not be null");
    set_mean(mean);
    set_stddev(stddev);
    double loglike = this->evaluate().loglike;
    if (!std::isfinite(loglike)) {
        throw std::invalid_argument(this->str() + " has non-finite loglike=" + std::to_string(loglike)
                                    + " on init");
    }
}

GaussianPrior::~GaussianPrior(){};

PriorEvaluation GaussianPrior::evaluate(bool calc_jacobians, bool normalize) const {
    double residual
            = ((_transformed ? _param->get_value_transformed() : _param->get_value()) - _mean) / _stddev;
    if (!std::isfinite(residual)) {
        throw std::runtime_error(this->str()
                                 + ".evaluate() got non-finite residual=" + std::to_string(residual));
    }
    double prior = -residual * residual / 2.;
    if (normalize) prior += this->get_loglike_const_terms()[0];
    std::map<ParamBaseCRef, std::vector<double>> jacobians = {};
    // dresidual/dx = d((x - mean)/stddev)/dx = 1/stddev
    if (calc_jacobians) {
        jacobians[*_param] = {1.0 / (_stddev * (_transformed ? 1.0 : _param->get_transform_derivative()))};
    }

    return PriorEvaluation{prior, {residual}, jacobians};
}

std::vector<double> GaussianPrior::get_loglike_const_terms() const {
    return {LOG_1_M_LOG_SQRT_2_PI - log(_stddev)};
}

double GaussianPrior::get_mean() const { return _mean; };
const ParamBase& GaussianPrior::get_param() const { return *_param; };
double GaussianPrior::get_stddev() const { return _stddev; };
bool GaussianPrior::get_transformed() const { return _transformed; };

void GaussianPrior::set_mean(double mean) {
    if (!std::isfinite(mean))
        throw std::invalid_argument("mean=" + to_string_float(mean) + " must be finite");
    _mean = mean;
}
void GaussianPrior::set_stddev(double stddev) {
    if (!(std::isfinite(stddev) && (stddev > 0))) {
        throw std::invalid_argument("stddev=" + to_string_float(stddev) + " must be finite and >0");
    }
    _stddev = stddev;
}
void GaussianPrior::set_transformed(bool transformed) { _transformed = transformed; };

size_t GaussianPrior::size() const { return 1; };

std::string GaussianPrior::repr(bool name_keywords, std::string_view namespace_separator) const {
    return type_name_str<GaussianPrior>(false, namespace_separator) + "(" + (name_keywords ? "param=" : "")
           + _param->repr(name_keywords, namespace_separator) + ", " + (name_keywords ? "mean=" : "")
           + to_string_float(_mean) + ", " + (name_keywords ? "stddev=" : "") + to_string_float(_stddev)
           + ", " + (name_keywords ? "transformed=" : "") + std::to_string(_transformed) + ")";
}

std::string GaussianPrior::str() const {
    return type_name_str<GaussianPrior>(true) + "(param=" + _param->str() + ", mean=" + to_string_float(_mean)
           + ", stddev=" + std::to_string(_stddev) + ", transformed=" + to_string_float(_transformed) + ")";
}
}  // namespace lsst::gauss2d::fit
