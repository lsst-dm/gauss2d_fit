#ifndef GAUSS2D_FIT_GSLSERSICMIXINTERPOLATOR_H
#define GAUSS2D_FIT_GSLSERSICMIXINTERPOLATOR_H

#ifdef GAUSS2D_FIT_HAS_GSL

#include <memory>

#include "gsl.h"
#include "gslinterpolator.h"
#include "sersicmix.h"

#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>

namespace gauss2d::fit {

/**
 * A SersicMixInterpolator that uses GSL interpolators between knots.
 */
class GSLSersicMixInterpolator : public SersicMixInterpolator {
private:
    mutable double _final_correction = 1.;
    const unsigned short _order;
    std::vector<std::pair<std::unique_ptr<GSLInterpolator>, std::unique_ptr<GSLInterpolator>>> _interps;

public:
    bool correct_final_integral = true;
    const GSLInterpType interp_type;
    /// The knot positions and values.
    const std::vector<SersicMixValues>& knots;

    double get_final_correction() const;

    std::vector<IntegralSize> get_integralsizes(double sersicindex) const override;
    std::vector<IntegralSize> get_integralsizes_derivs(double sersicindex) const override;

    unsigned short get_order() const override;

    const double sersicindex_min;
    const double sersicindex_max;

    std::string repr(bool name_keywords = false) const override;
    std::string str() const override;

    explicit GSLSersicMixInterpolator(
        unsigned short order = SERSICMIX_ORDER_DEFAULT,
        const GSLInterpType interp_type=GSLInterpType::cspline);
    ~GSLSersicMixInterpolator();
};

}  // namespace gauss2d::fit

#endif // GAUSS2D_FIT_HAS_GSL
#endif // GAUSS2D_FIT_GSLSERSICMIXINTERPOLATOR_H