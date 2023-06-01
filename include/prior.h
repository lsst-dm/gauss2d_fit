#ifndef GAUSS2D_FIT_PRIOR_H
#define GAUSS2D_FIT_PRIOR_H

#include <map>

#include "gauss2d/object.h"

#include "math.h"
#include "param_defs.h"

namespace gauss2d::fit {

/**
 * Results from the evaluation of a prior probability function.
 */
struct PriorEvaluation {
public:
    double loglike;
    std::map<ParamBaseCRef, std::vector<double>> jacobians = {};
    std::vector<double> residuals = {};
};

/**
 * Interface for a prior probability function.
 */
class Prior : public Object {
public:
    /// Return the constant terms of the log likelihood (dependent on stddevs only)
    virtual std::vector<double> get_loglike_const_terms() const = 0;
    /**
     * Evaluate the log likelihood and residual-dependent terms.
     *
     * @param calc_jacobians Whether to compute the Jacobian and residual terms.
     * @param normalize_loglike Whether to add the constant (sigma-dependent) factors to the log likehood.
     * @return The result of the evaluation.
     */
    virtual PriorEvaluation evaluate(bool calc_jacobians = false, bool normalize_loglike = false) const = 0;

    virtual size_t size() const = 0;

    virtual ~Prior() = default;
};
}  // namespace gauss2d::fit

#endif
