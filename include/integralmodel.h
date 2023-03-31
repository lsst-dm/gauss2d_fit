#ifndef GAUSS2D_FIT_INTEGRALMODEL_H
#define GAUSS2D_FIT_INTEGRALMODEL_H

#include <memory>
#include <set>

#include "channel.h"
#include "parametric.h"
#include "parametricmodel.h"

namespace gauss2d::fit
{

/*
    An IntegralModel describes the total integrated surface density of a
    2D distribution. It is intended to define a Component but can
    be used otherwise.

    It is generally not intended for IntegralModels to be shared amongst
    Components or otherwise, but this is not guaranteed to be enforced.
*/
class IntegralModel : public Parametric
{
public:
    virtual std::set<std::reference_wrapper<const Channel>> get_channels() const = 0;
    virtual double get_integral(const Channel & channel) const = 0;
    // Should be 1 for a linear model and ?? for non-linear
    // Needed to re-compute derivative factors when fitting non-linear models
    virtual std::vector<std::pair<ParamBaseCRef, extra_param_factor_values>> get_integral_derivative_factors(
        const Channel & channel) const = 0;
};

inline bool operator < ( const IntegralModel &lhs, const IntegralModel &rhs) { return &lhs < &rhs; }
//const bool operator == ( const IntegralModel &m ) const { return &(*this) == &m; };
//const bool operator != ( const IntegralModel &m ) const { return &(*this) != &m; };


} // namespace gauss2d::fit

#endif
