#ifndef GAUSS2D_FIT_INTEGRALMODEL_H
#define GAUSS2D_FIT_INTEGRALMODEL_H

#include <memory>
#include <set>

#include "channel.h"
#include "parametric.h"
#include "parametricmodel.h"

namespace gauss2d::fit
{
/**
 * @brief A Parametric model for the integral of a 2D distribution.
 *
 * IntegralModel implementations can use as many Parameter classes as needed.
 * The model should be efficiently differentiable in the first order, via
 * finite differencing if analytical evaluation is not possible.
 */
class IntegralModel : public Parametric
{
public:
    /// Get the channels this instance is applicable for
    virtual std::set<std::reference_wrapper<const Channel>> get_channels() const = 0;
    /// Get the value of the integral in a single Channel
    virtual double get_integral(const Channel & channel) const = 0;
    /**
     * @brief Return the partial derivative of the model w.r.t. each metaparameter.
     *
     * This is generally needed only for nonlinear Parameter instances.
     *
     * @param channel The Channel to return derivatives for.
     * @return A vector of Parameter/derivative pairs, with values ordered as
     *         specified in GaussianEvaluator (L, sigma_x, sigma_y).
     */
    virtual std::vector<std::pair<ParamBaseCRef, ExtraParamFactorValues>> get_integral_derivative_factors(
        const Channel & channel) const = 0;
};

inline bool operator < ( const IntegralModel &lhs, const IntegralModel &rhs) { return &lhs < &rhs; }
// TODO: These aren't necessary, but should they be included?
//const bool operator == ( const IntegralModel &m ) const { return &(*this) == &m; };
//const bool operator != ( const IntegralModel &m ) const { return &(*this) != &m; };


} // namespace gauss2d::fit

#endif
