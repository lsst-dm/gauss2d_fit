#ifndef GAUSS2D_FIT_SERSICMIX_H
#define GAUSS2D_FIT_SERSICMIX_H

#include <array>
#include <stdexcept>
#include <vector>

#include "gauss2d/object.h"

namespace gauss2d
{
namespace fit
{

static const short SERSICMIX_ORDER_DEFAULT = 4;

class IntegralSize : public Object {
public:
    const double integral;
    const double sigma;

    std::string str() const override;
    IntegralSize(double integral_, double sigma_);
};

class SersicMixInterpolator : public Object {
public:
    virtual std::vector<IntegralSize> get_integralsizes(double sersicindex) const = 0;
    virtual std::vector<IntegralSize> get_integralsizes_derivs(double sersicindex) const = 0;

    virtual unsigned short get_order() const = 0;

    virtual ~SersicMixInterpolator() {};
};

class SersicMixValues : public Object {
public:
    const double sersicindex;
    const std::vector<IntegralSize> values;

    std::string str() const override;
    SersicMixValues(double sersicindex_, std::vector<IntegralSize> values_);
};

inline bool operator < ( const SersicMixValues &lhs, const SersicMixValues &rhs) {
    return lhs.sersicindex < rhs.sersicindex;
}
inline bool operator < ( const SersicMixValues &lhs, double x) { return lhs.sersicindex < x; }
inline bool operator < ( double x, const SersicMixValues &rhs) { return x < rhs.sersicindex; }

std::vector<SersicMixValues> get_sersic_mix_knots_copy(unsigned short order);

const std::vector<SersicMixValues> & get_sersic_mix_knots(unsigned short order);

const std::vector<SersicMixValues> & get_sersic_mix_knots_order4();
const std::vector<SersicMixValues> & get_sersic_mix_knots_order8();

} // namespace fit
} // namespace gauss2d

#endif
