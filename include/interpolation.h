#ifndef GAUSS2DFIT_INTERPOLATION_H
#define GAUSS2DFIT_INTERPOLATION_H

namespace gauss2d::fit {

enum class InterpType {
    linear,      ///< Linear interpolation.
    polynomial,  ///< Polynomial interpolation.
    cspline,     ///< Cubic spline interpolation.
    akima,       ///< Akima spline with natural boundary conditions.
};

}  // namespace gauss2d::fit

#endif  // GAUSS2DFIT_INTERPOLATION_H
