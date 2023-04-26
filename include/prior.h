#ifndef GAUSS2D_FIT_PRIOR_H
#define GAUSS2D_FIT_PRIOR_H

#include <map>

#include "gauss2d/object.h"

#include "param_defs.h"

const double LOG_1 = log(1);
const double SQRT_2_PI = sqrt(2.*M_PI);

namespace gauss2d::fit {
/**
 * Results from the evaluation of a prior probability function.
 */
struct PriorEvaluation {
public:
    double loglike;
    std::map<ParamBaseCRef, double> jacobians = {};
    std::vector<double> residuals = {};
};

/**
 * Interface for a prior probability function.
 */
class Prior : public Object {
public:
    virtual double get_loglike_const_term() const = 0;
    virtual PriorEvaluation evaluate(bool calc_jacobians=false, bool normalize_loglike=false) const = 0;

    virtual size_t size() const = 0;

    virtual ~Prior() = default;
};
}  // namespace gauss2d::fit

#endif
