#ifndef LSST_GAUSS2D_FIT_PRIOR_H
#define LSST_GAUSS2D_FIT_PRIOR_H

#include <map>

#include "lsst/gauss2d/object.h"

#include "math.h"
#include "param_defs.h"

namespace lsst::gauss2d::fit{

/**
 * Results from the evaluation of a prior probability function.
 */
class PriorEvaluation : public Object {
public:
    typedef std::map<ParamBaseCRef, std::vector<double>> Jacobians;

    double loglike;
    std::vector<double> residuals;
    Jacobians jacobians;

    double compute_dloglike_dx(const ParamBase& param, bool transformed = true) const;

    std::string repr(bool name_keywords = false,  std::string_view namespace_separator = Object::CC_NAMESPACE_SEPARATOR) const override;
    std::string str() const override;

    PriorEvaluation(double loglike, std::vector<double> residuals = {}, Jacobians jacobians = {},
                    bool check_size = true);
    ~PriorEvaluation();
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
    virtual PriorEvaluation evaluate(bool calc_jacobians = false, bool normalize_loglike = true) const = 0;

    virtual size_t size() const = 0;

    virtual ~Prior() = default;
};
}  // namespace lsst::gauss2d::fit

#endif
