#ifndef GAUSS2D_FIT_PARAMETRICMODEL_H
#define GAUSS2D_FIT_PARAMETRICMODEL_H

#include "gauss2d/evaluate.h"
#include "gauss2d/gaussian.h"

#include "channel.h"
#include "parametric.h"

namespace gauss2d
{
namespace fit
{
typedef std::vector<std::array<size_t, gauss2d::N_EXTRA_MAP>> extra_param_map;
typedef std::vector<std::array<double, gauss2d::N_EXTRA_FACTOR>> extra_param_factors;
typedef std::vector<std::array<size_t, gauss2d::N_PARAMS>> grad_param_map;
typedef std::vector<std::array<double, gauss2d::N_PARAMS>> grad_param_factors;

typedef std::map<ParamBaseCRef, size_t> ParameterMap;

class ParametricModel : public Parametric
{
public:
    virtual void add_extra_param_map(const Channel & channel, extra_param_map & map, ParameterMap & offsets
        ) const = 0;
    virtual void add_extra_param_factors(const Channel & channel, extra_param_factors & factors) const = 0;
    virtual void add_grad_param_map(const Channel & channel, grad_param_map & map, ParameterMap & offsets
        ) const = 0;
    virtual void add_grad_param_factors(const Channel & channel, grad_param_factors & factor) const = 0;
    
    virtual std::unique_ptr<const gauss2d::Gaussians> get_gaussians(const Channel & channel) const = 0;
};

} // namespace fit
} // namespace gauss2d

#endif