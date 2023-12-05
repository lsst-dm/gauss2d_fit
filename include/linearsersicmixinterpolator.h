#ifndef GAUSS2D_FIT_LINEARSERSICMIXINTERPOLATOR_H
#define GAUSS2D_FIT_LINEARSERSICMIXINTERPOLATOR_H

#include "sersicmix.h"

namespace gauss2d::fit {

/**
 * A SersicMixInterpolator that uses linear interpolation between knots.
 */
class LinearSersicMixInterpolator : public SersicMixInterpolator {
private:
    const unsigned short _order;

public:
    /// The knot positions and values.
    const std::vector<SersicMixValues>& knots;

    std::vector<IntegralSize> get_integralsizes(double sersicindex) const override;
    std::vector<IntegralSize> get_integralsizes_derivs(double sersicindex) const override;

    InterpType get_interptype() const override;
    unsigned short get_order() const override;

    const double sersicindex_min;
    const double sersicindex_max;

    std::string repr(bool name_keywords = false) const override;
    std::string str() const override;

    explicit LinearSersicMixInterpolator(unsigned short order = SERSICMIX_ORDER_DEFAULT);
    ~LinearSersicMixInterpolator();
};

}  // namespace gauss2d::fit

#endif
