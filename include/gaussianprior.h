#ifndef GAUSS2D_FIT_LSQPRIOR_H
#define GAUSS2D_FIT_LSQPRIOR_H

#include <memory>

#include "gauss2d/object.h"

#include "param_defs.h"
#include "prior.h"

namespace gauss2d::fit {
/**
 * A 1D Gaussian prior for a single Parameter.
 */
class GaussianPrior : public Prior {
private:
    std::shared_ptr<const ParamBase> _param;
    double _mean;
    double _stddev;
    bool _transformed;

public:
    /**
     * Construct a GaussianPrior from a Parameter and mean/std. deviation.
     *
     * @param param The ParamBase to compute a prior for.
     * @param mean The mean value of the prior.
     * @param stddev The standard deviation of the prior.
     * @param transformed Whether the prior is based on the transformed value of param.
     */
    GaussianPrior(std::shared_ptr<const ParamBase> param, double mean, double stddev, bool transformed);
    /**
     * Evaluate the log likelihood and residual-dependent terms.
     *
     * @param calc_jacobians Whether to compute the Jacobian and residual terms.
     * @param normalize_loglike Whether to add the constant (sigma-dependent) factors to the log likehood.
     * @return The result of the evaluation.
     */
    PriorEvaluation evaluate(bool calc_jacobians=false, bool normalize_loglike=false) const;

    const ParamBase & get_param() const;

    double get_loglike_const_term() const override;
    double get_mean() const;
    double get_stddev() const;
    bool get_transformed() const;

    void set_mean(double mean);
    void set_stddev(double stddev);
    void set_transformed(bool transformed);

    size_t size() const override;

    std::string repr(bool name_keywords = false) const override;
    std::string str() const override;

    ~GaussianPrior();
};
}  // namespace gauss2d::fit

#endif
