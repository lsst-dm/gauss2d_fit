#ifndef LSST_GAUSS2D_FIT_GSLSERSICMIXINTERPOLATOR_H
#define LSST_GAUSS2D_FIT_GSLSERSICMIXINTERPOLATOR_H

#ifdef LSST_GAUSS2D_FIT_HAS_GSL

#include <memory>

#include "gsl.h"
#include "gslinterpolator.h"
#include "sersicmix.h"

#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>

namespace lsst::gauss2d::fit {

/**
 * A SersicMixInterpolator that uses GSL interpolators between knots.
 */
class GSLSersicMixInterpolator : public SersicMixInterpolator {
private:
    mutable double _final_correction = 1.;
    const InterpType _interp_type;
    const unsigned short _order;
    std::vector<std::pair<std::unique_ptr<GSLInterpolator>, std::unique_ptr<GSLInterpolator>>> _interps;

public:
    bool correct_final_integral = true;
    static constexpr InterpType INTERPTYPE_DEFAULT = InterpType::cspline;
    /// The knot positions and values.
    const std::vector<SersicMixValues>& knots;

    /// Get the multiplicative factor required to adjust the integral for the final order
    /// component such that the sum of all integral factors is unity (normalized).
    double get_final_correction() const;

    std::vector<IntegralSize> get_integralsizes(double sersicindex) const override;
    std::vector<IntegralSize> get_integralsizes_derivs(double sersicindex) const override;

    InterpType get_interptype() const override;
    unsigned short get_order() const override;

    const double sersicindex_min;
    const double sersicindex_max;

    std::string repr(bool name_keywords = false,
                     std::string_view namespace_separator = Object::CC_NAMESPACE_SEPARATOR) const override;
    std::string str() const override;

    explicit GSLSersicMixInterpolator(unsigned short order = SERSICMIX_ORDER_DEFAULT,
                                      InterpType interp_type = INTERPTYPE_DEFAULT);
    ~GSLSersicMixInterpolator();
};

}  // namespace lsst::gauss2d::fit

#endif  // LSST_GAUSS2D_FIT_HAS_GSL
#endif  // LSST_GAUSS2D_FIT_GSLSERSICMIXINTERPOLATOR_H