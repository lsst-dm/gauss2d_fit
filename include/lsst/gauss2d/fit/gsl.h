#ifndef LSST_GAUSS2D_FIT_GSL_H
#define LSST_GAUSS2D_FIT_GSL_H

#ifdef LSST_GAUSS2D_FIT_HAS_GSL

#include <unordered_map>

#include <gsl/gsl_interp.h>

#include "interpolation.h"

namespace lsst::gauss2d::fit{

/**
 * See GSL docs for 1D interpolation types, or (as of 2.7):
 * https://www.gnu.org/software/gsl/doc/html/interp.html#d-interpolation-types
 */
static const std::unordered_map<InterpType, const gsl_interp_type*> GSLInterpTypes{
        {InterpType::linear, gsl_interp_linear},
        {InterpType::polynomial, gsl_interp_polynomial},
        {InterpType::cspline, gsl_interp_cspline},
        {InterpType::akima, gsl_interp_akima}};
}  // namespace lsst::gauss2d::fit

#endif  // LSST_GAUSS2D_FIT_HAS_GSL
#endif  // GAUSS2D_FIT_GSL_H
