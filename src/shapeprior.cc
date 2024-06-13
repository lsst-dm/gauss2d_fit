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

#include "lsst/gauss2d/ellipse.h"
#include "lsst/gauss2d/to_string.h"
#include "lsst/gauss2d/type_name.h"

#include "lsst/gauss2d/fit/param_defs.h"
#include "lsst/gauss2d/fit/shapeprior.h"

namespace lsst::gauss2d::fit {
ShapePriorOptions::ShapePriorOptions(double delta_jacobian, double size_maj_floor, double axrat_floor)
        : _delta_jacobian(delta_jacobian), _size_maj_floor(size_maj_floor), _axrat_floor(axrat_floor) {
    this->check_delta_jacobian(_delta_jacobian, true);
    this->check_size_maj_floor(_size_maj_floor, true);
    this->check_axrat_floor(_axrat_floor, true);
}

bool ShapePriorOptions::check_delta_jacobian(double delta_jacobian, bool do_throw) {
    if (!std::isfinite(delta_jacobian)) {
        if (do_throw) {
            throw std::invalid_argument("delta_jacobian=" + std::to_string(delta_jacobian)
                                        + " must be finite");
        }
        return false;
    }
    return true;
}

bool ShapePriorOptions::check_size_maj_floor(double size_maj_floor, bool do_throw) {
    if (!(std::isfinite(size_maj_floor) && (size_maj_floor > 0))) {
        if (do_throw) {
            throw std::invalid_argument("size_maj_floor=" + std::to_string(size_maj_floor)
                                        + " must be >0 and finite");
        }
        return false;
    }
    return true;
}

bool ShapePriorOptions::check_axrat_floor(double axrat_floor, bool do_throw) {
    if (!(std::isfinite(axrat_floor) && (axrat_floor > 0))) {
        if (do_throw) {
            throw std::invalid_argument("axrat_floor=" + std::to_string(axrat_floor)
                                        + " must be >0 and finite");
        }
        return false;
    }
    return true;
}

double ShapePriorOptions::get_delta_jacobian() const { return _delta_jacobian; }

double ShapePriorOptions::get_size_maj_floor() const { return _size_maj_floor; }

double ShapePriorOptions::get_axrat_floor() const { return _axrat_floor; }

void ShapePriorOptions::set_delta_jacobian(double delta_jacobian) {
    check_delta_jacobian(delta_jacobian, true);
    _delta_jacobian = delta_jacobian;
}

void ShapePriorOptions::set_size_maj_floor(double size_maj_floor) {
    check_size_maj_floor(size_maj_floor, true);
    this->_size_maj_floor = size_maj_floor;
}

void ShapePriorOptions::set_axrat_floor(double axrat_floor) {
    check_axrat_floor(axrat_floor, true);
    this->_axrat_floor = axrat_floor;
}

std::string ShapePriorOptions::repr(bool name_keywords, std::string_view namespace_separator) const {
    return type_name_str<ShapePriorOptions>(false, namespace_separator) + "("
           + (name_keywords ? "delta_jacobian=" : "") + to_string_float(_delta_jacobian) + ", "
           + (name_keywords ? "size_maj_floor=" : "") + to_string_float(_size_maj_floor) + ", "
           + (name_keywords ? "axrat_floor=" : "") + to_string_float(_axrat_floor) + ")";
}

std::string ShapePriorOptions::str() const {
    return type_name_str<ShapePriorOptions>(true) + "(delta_jacobian=" + to_string_float(_delta_jacobian)
           + ", size_maj_floor=" + to_string_float(_size_maj_floor)
           + ", axrat_floor=" + to_string_float(_axrat_floor) + ")";
}

ShapePrior::ShapePrior(std::shared_ptr<const ParametricEllipse> ellipse,
                       std::shared_ptr<ParametricGaussian1D> prior_size,
                       std::shared_ptr<ParametricGaussian1D> prior_axrat,
                       std::shared_ptr<ShapePriorOptions> options)
        : _ellipse(std::move(ellipse)),
          _prior_size(std::move(prior_size)),
          _prior_axrat(std::move(prior_axrat)),
          _options(std::move((options == nullptr) ? std::make_shared<ShapePriorOptions>() : options)) {
    if (_ellipse == nullptr) throw std::invalid_argument("ellipse param must not be null");
    double loglike = this->evaluate().loglike;
    if (!std::isfinite(loglike)) {
        throw std::invalid_argument(this->str() + " has non-finite loglike=" + std::to_string(loglike)
                                    + " on init");
    }
}

ShapePrior::~ShapePrior() {}

PriorEvaluation ShapePrior::evaluate(bool calc_jacobians, bool normalize) const {
    const auto ellipse = lsst::gauss2d::Ellipse(this->_ellipse->get_size_x(), this->_ellipse->get_size_y(),
                                                this->_ellipse->get_rho());
    const auto ellipse_major = lsst::gauss2d::EllipseMajor(ellipse);

    double size_maj = ellipse_major.get_r_major();
    double axrat = size_maj == 0 ? 0 : ellipse_major.get_axrat();

    if (!(axrat >= 0)) {
        throw std::runtime_error(this->str() + " got invalid axrat=" + std::to_string(axrat));
    }

    auto result = PriorEvaluation(0, {}, {}, false);

    if (this->_prior_size != nullptr) {
        double size_maj_floor = this->_options->get_size_maj_floor();
        size_maj = sqrt(size_maj * size_maj * axrat + size_maj_floor * size_maj_floor);
        size_maj = _prior_size->get_mean_parameter().get_transform().forward(size_maj);
        double sigma = _prior_size->get_stddev();
        double residual = (size_maj - _prior_size->get_mean_parameter().get_value_transformed()) / sigma;
        if (!std::isfinite(residual)) {
            throw std::runtime_error(this->str() + ".evaluate() got non-finite size residual="
                                     + std::to_string(residual) + " from size_maj=" + std::to_string(size_maj)
                                     + ", _prior_size=" + _prior_size->str());
        }
        result.residuals.emplace_back(residual);
        result.loglike += normalize ? logpdf_norm(residual, 1.0) : -residual * residual / 2.;
    }

    if (this->_prior_axrat != nullptr) {
        double axrat_floor = this->_options->get_axrat_floor();
        axrat = sqrt(axrat_floor * axrat_floor + axrat * axrat);
        if (axrat > 1) {
            axrat = 1;
        } else if (!(axrat >= 0)) {
            throw std::runtime_error(this->str() + " got invalid axrat=" + std::to_string(axrat)
                                     + " from axrat_floor=" + std::to_string(axrat_floor));
        }
        const auto& transform_param = _prior_axrat->get_mean_parameter().get_transform();
        axrat = transform_param.forward(axrat);
        double sigma = _prior_axrat->get_stddev();
        double residual = (axrat - _prior_axrat->get_mean_parameter().get_value_transformed()) / sigma;
        if (!std::isfinite(residual)) {
            auto str = this->str();
            throw std::runtime_error(str + ".evaluate() got non-finite axrat residual="
                                     + std::to_string(residual) + " from axrat=" + std::to_string(axrat)
                                     + ", _prior_axrat=" + _prior_axrat->str());
        }
        // RuntimeError('Infinite axis ratio prior residual from q={axrat} and mean, std, f
        // 'logit stretch divisor = {self.axrat_params}')
        result.residuals.emplace_back(residual);
        result.loglike += normalize ? logpdf_norm(residual, 1.0) : -residual * residual / 2.;
    }

    if (calc_jacobians) {
        ParamRefs params = nonconsecutive_unique(_ellipse->get_parameters_new());
        for (auto& paramref : params) {
            auto& param = paramref.get();
            double value_init = param.get_value();
            double value_trans = param.get_value_transformed();
            double delta = finite_difference_param(param, _options->get_delta_jacobian());
            double value_new = param.get_value_transformed();

            std::vector<double> residuals_old;
            std::vector<double> residuals_new;
            try {
                param.set_value_transformed(value_trans + delta / 2.);
                residuals_new = this->evaluate(false, normalize).residuals;
                param.set_value_transformed(value_trans - delta / 2.);
                residuals_old = this->evaluate(false, normalize).residuals;
            } catch (std::exception& e) {
                param.set_value_transformed(value_new);
                residuals_new = this->evaluate(false, normalize).residuals;
                residuals_old = result.residuals;
            }
            // Return to the old value
            param.set_value(value_init);
            value_new = param.get_value();
            if (value_new != value_init) {
                throw std::runtime_error(this->str() + " could not return " + param.str()
                                         + " to original value=" + to_string_float(value_init) + " (stuck at "
                                         + to_string_float(value_new) + "); check limits");
            }
            std::vector<double> jacobians(result.residuals.size());
            for (size_t idx = 0; idx < result.residuals.size(); ++idx) {
                jacobians[idx] = (residuals_new[idx] - residuals_old[idx]) / delta;
            }
            result.jacobians[paramref] = jacobians;
        }
    }

    return result;
}

std::shared_ptr<ParametricGaussian1D> ShapePrior::get_prior_size() const { return _prior_size; }

std::shared_ptr<ParametricGaussian1D> ShapePrior::get_prior_axrat() const { return _prior_axrat; }

std::vector<double> ShapePrior::get_loglike_const_terms() const {
    return {this->_prior_size == nullptr ? 0 : (LOG_1 - log(this->_prior_size->get_stddev() * SQRT_2_PI)),
            this->_prior_axrat == nullptr ? 0 : (LOG_1 - log(this->_prior_axrat->get_stddev() * SQRT_2_PI))};
}

void ShapePrior::set_prior_size(std::shared_ptr<ParametricGaussian1D> prior_size) {
    _prior_size = prior_size;
}

void ShapePrior::set_prior_axrat(std::shared_ptr<ParametricGaussian1D> prior_axrat) {
    _prior_axrat = prior_axrat;
}

size_t ShapePrior::size() const { return (this->_prior_size != nullptr) + (this->_prior_axrat != nullptr); };

std::string ShapePrior::repr(bool name_keywords, std::string_view namespace_separator) const {
    return type_name_str<ShapePrior>(false, namespace_separator) + "(" + (name_keywords ? "ellipse=" : "")
           + _ellipse->repr(name_keywords, namespace_separator) + ", " + (name_keywords ? "prior_size=" : "")
           + repr_ptr(_prior_size.get(), name_keywords, namespace_separator) + ", "
           + (name_keywords ? "prior_axrat=" : "")
           + repr_ptr(_prior_axrat.get(), name_keywords, namespace_separator) + ", "
           + (name_keywords ? "options=" : "") + _options->repr(name_keywords, namespace_separator) + ")";
}

std::string ShapePrior::str() const {
    return type_name_str<ShapePrior>(true) + "(ellipse=" + _ellipse->str()
           + ", prior_size=" + str_ptr(_prior_size.get()) + ", prior_axrat=" + str_ptr(_prior_axrat.get())
           + ", options=" + _options->str() + ")";
}

}  // namespace lsst::gauss2d::fit
