#ifndef LSST_GAUSS2D_FIT_LINEARSERSICMIXINTERPOLATOR_H
#define LSST_GAUSS2D_FIT_LINEARSERSICMIXINTERPOLATOR_H

#include "sersicmix.h"

namespace lsst::gauss2d::fit {

/**
 * A SersicMixInterpolator that uses linear interpolation between knots.
 */
class LinearSersicMixInterpolator : public SersicMixInterpolator {
public:
    explicit LinearSersicMixInterpolator(unsigned short order = SERSICMIX_ORDER_DEFAULT);
    ~LinearSersicMixInterpolator();

    std::vector<IntegralSize> get_integralsizes(double sersicindex) const override;
    std::vector<IntegralSize> get_integralsizes_derivs(double sersicindex) const override;
    /// The knot positions and values.
    const std::vector<SersicMixValues>& get_knots() const;

    InterpType get_interptype() const override;
    unsigned short get_order() const override;
    double get_sersicindex_min() const;
    double get_sersicindex_max() const;

    std::string repr(bool name_keywords = false,
                     std::string_view namespace_separator = Object::CC_NAMESPACE_SEPARATOR) const override;
    std::string str() const override;

private:
    const unsigned short _order;
    const std::vector<SersicMixValues>& _knots;
    const double _sersicindex_min;
    const double _sersicindex_max;
};

}  // namespace lsst::gauss2d::fit

#endif
