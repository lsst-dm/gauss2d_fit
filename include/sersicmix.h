#ifndef GAUSS2D_FIT_SERSICMIX_H
#define GAUSS2D_FIT_SERSICMIX_H

#include <array>
#include <memory>
#include <stdexcept>
#include <vector>

#include "gauss2d/object.h"

namespace gauss2d::fit {

static const unsigned short SERSICMIX_ORDER_DEFAULT = 4;

/**
 * A pair of integral - size values for a Gaussian (sub)Component.
 */
class IntegralSize : public Object {
public:
    const double integral;
    const double sigma;

    std::string repr(bool name_keywords = false) const override;
    std::string str() const override;
    IntegralSize(double integral_, double sigma_);
};

/**
 * An interpolator that returns IntegralSize vectors for a given Sersic index.
 */
class SersicMixInterpolator : public Object {
public:
    /**
     * Get the vector of IntegralSize values for a given Sersic index.
     *
     * @param sersicindex The Sersic index value.
     * @return The vector of IntegralSize values for sersicindex.
     */
    virtual std::vector<IntegralSize> get_integralsizes(double sersicindex) const = 0;
    virtual std::vector<IntegralSize> get_integralsizes_derivs(double sersicindex) const = 0;

    virtual unsigned short get_order() const = 0;
};

const std::shared_ptr<const SersicMixInterpolator> get_sersic_mix_interpolator_default(
        unsigned short order=SERSICMIX_ORDER_DEFAULT);

/**
 * A vector of IntegralSize values for a given Sersic index.
 */
class SersicMixValues : public Object {
public:
    const double sersicindex;
    const std::vector<IntegralSize> values;

    std::string repr(bool name_keywords = false) const override;
    std::string str() const override;
    /**
     * Construct a SersicMixValues.
     *
     * @param sersicindex The value of the Sersic index.
     * @param values The vector of IntegralSize values for the Sersic index.
     */
    SersicMixValues(double sersicindex, std::vector<IntegralSize> values);
};

inline bool operator<(const SersicMixValues &lhs, const SersicMixValues &rhs) {
    return lhs.sersicindex < rhs.sersicindex;
}
inline bool operator<(const SersicMixValues &lhs, double x) { return lhs.sersicindex < x; }
inline bool operator<(double x, const SersicMixValues &rhs) { return x < rhs.sersicindex; }

std::vector<SersicMixValues> get_sersic_mix_knots_copy(unsigned short order);

const std::vector<SersicMixValues> &get_sersic_mix_knots(unsigned short order);

const std::vector<SersicMixValues> &get_sersic_mix_knots_order4();
const std::vector<SersicMixValues> &get_sersic_mix_knots_order8();

}  // namespace gauss2d::fit

#endif
