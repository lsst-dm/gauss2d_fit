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

#ifdef GAUSS2D_FIT_HAS_GSL

#ifndef GAUSS2D_FIT_GSLINTERPOLATOR_H
#define GAUSS2D_FIT_GSLINTERPOLATOR_H

#include "gsl.h"
#include "interpolation.h"

#include "gauss2d/object.h"

#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>

namespace gauss2d::fit {

/**
 * A wrapper for GSL 1D interpolators.
 */
class GSLInterpolator : public Object {
private:
    gsl_interp_accel* _acc;
    gsl_spline* _spline;
    size_t _n_knots;
    std::vector<double> _x;
    std::vector<double> _y;

public:
    const InterpType interp_type;
    static constexpr InterpType INTERPTYPE_DEFAULT = InterpType::cspline;

    /// Get the interpolant value for a knot of the given index
    double get_knot_x(size_t idx) const;
    /// Get the interpolated function value for a knot of the given index
    double get_knot_y(size_t idx) const;

    /// Get the interpolated function value for a given interpolant value
    double eval(double x) const;
    /// Get the derivative of the interpolated function value for a given interpolant value
    double eval_deriv(double x) const;
    /// Get the number of knots
    size_t size() const;

    std::string repr(bool name_keywords = false) const override;
    std::string str() const override;

    explicit GSLInterpolator(std::vector<double> x, std::vector<double> y,
                             InterpType interp_type = INTERPTYPE_DEFAULT);
    ~GSLInterpolator();
};

}  // namespace gauss2d::fit

#endif  // GAUSS2D_FIT_GSLINTERPOLATOR_H
#endif  // GAUSS2D_FIT_HAS_GSL
