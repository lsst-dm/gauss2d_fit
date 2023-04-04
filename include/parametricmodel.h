#ifndef GAUSS2D_FIT_PARAMETRICMODEL_H
#define GAUSS2D_FIT_PARAMETRICMODEL_H

#include "gauss2d/evaluate.h"
#include "gauss2d/gaussian.h"

#include "channel.h"
#include "parametric.h"

namespace gauss2d::fit
{
typedef std::vector<std::array<size_t, gauss2d::N_EXTRA_MAP>> ExtraParamMap;
typedef std::array<double, gauss2d::N_EXTRA_FACTOR> ExtraParamFactorValues;
typedef std::vector<ExtraParamFactorValues> ExtraParamFactors;
typedef std::vector<std::array<size_t, gauss2d::N_PARAMS_GAUSS2D>> GradParamMap;
typedef std::vector<std::array<double, gauss2d::N_PARAMS_GAUSS2D>> GradParamFactors;

typedef std::map<ParamBaseCRef, size_t> ParameterMap;

class ParametricModel : public Parametric
{
public:
    virtual void add_extra_param_map(
            const Channel & channel, ExtraParamMap & map_extra, const GradParamMap & map_grad, ParameterMap & offsets
    ) const = 0;
    virtual void add_extra_param_factors(const Channel & channel, ExtraParamFactors & factors) const = 0;
    virtual void add_grad_param_map(const Channel & channel, GradParamMap & map, ParameterMap & offsets
        ) const = 0;
    virtual void add_grad_param_factors(const Channel & channel, GradParamFactors & factors) const = 0;
    virtual void set_extra_param_factors(const Channel & channel, ExtraParamFactors & factors, size_t index) const = 0;
    virtual void set_grad_param_factors(const Channel & channel, GradParamFactors & factors, size_t index) const = 0;

    virtual std::unique_ptr<const gauss2d::Gaussians> get_gaussians(const Channel & channel) const = 0;
    virtual size_t get_n_gaussians(const Channel & channel) const = 0;
};

} // namespace gauss2d::fit

#endif
