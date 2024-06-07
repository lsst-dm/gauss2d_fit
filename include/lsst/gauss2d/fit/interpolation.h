#ifndef LSST_GAUSS2D_FIT_INTERPOLATION_H
#define LSST_GAUSS2D_FIT_INTERPOLATION_H

namespace lsst::gauss2d::fit{

enum class InterpType {
    linear,      ///< Linear interpolation.
    polynomial,  ///< Polynomial interpolation.
    cspline,     ///< Cubic spline interpolation.
    akima,       ///< Akima spline with natural boundary conditions.
};

}  // namespace lsst::gauss2d::fit

#endif  // GAUSS2D_FIT_INTERPOLATION_H
