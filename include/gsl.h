#ifndef GAUSS2DFIT_GSL_H
#define GAUSS2DFIT_GSL_H

#ifdef GAUSS2D_FIT_HAS_GSL

#include <unordered_map>

#include <gsl/gsl_interp.h>

namespace gauss2d::fit {

/**
 * See GSL docs for 1D interpolation types, or (as of 2.7):
 * https://www.gnu.org/software/gsl/doc/html/interp.html#d-interpolation-types
 */
enum class GSLInterpType {
    linear,          ///< Linear interpolation.
    polynomial,      ///< Polynomial interpolation.
    cspline,         ///< Cubic spline with natural boundary conditions.
    akima,           ///< Non-rounded Akima spline with natural boundary conditions.
};

static const std::unordered_map<GSLInterpType, const gsl_interp_type*> GSLInterpTypes {
    {GSLInterpType::linear, gsl_interp_linear},
    {GSLInterpType::polynomial, gsl_interp_polynomial},
    {GSLInterpType::cspline, gsl_interp_cspline},
    {GSLInterpType::akima, gsl_interp_akima}
};
}  // namespace gauss2d::fit

#endif  // GAUSS2D_FIT_HAS_GSL
#endif  // GAUSS2DFIT_GSL_H
