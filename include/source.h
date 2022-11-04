#ifndef GAUSS2D_FIT_SOURCE_H
#define GAUSS2D_FIT_SOURCE_H

#include <memory>

#include "component.h"
#include "componentmixture.h"
#include "param_filter.h"

namespace gauss2d
{
namespace fit
{

/*
    A Source is a loose association of Components. The general intent is to
    group Components with the same (or slightly offset) centroids, but users
    can choose other systems if needed.

    Sources are also not meant to share Component pointers. Unfortunately,
    this can't be enforced in pybind11 bindings, so the constructor allows
    shared pointers, and there is no registry or other mechanism to 
    enforce a 1:1 Source:Component relationship. Caveat emptor.
*/
class Source : public ComponentMixture
{
private:
    Components _components = {};

public:
    void add_extra_param_map(const Channel & channel, extra_param_map & map_extra,
        const grad_param_map & map_grad, ParameterMap & offsets
    ) const override;
    void add_extra_param_factors(const Channel & channel, extra_param_factors & factors) const override;
    void add_grad_param_map(const Channel & channel, grad_param_map & map, ParameterMap & offsets
        ) const override;
    void add_grad_param_factors(const Channel & channel, grad_param_factors & factors) const override;
    
    Components get_components() const override;
    std::unique_ptr<const gauss2d::Gaussians> get_gaussians(const Channel & channel) const override;
    size_t get_n_gaussians(const Channel & channel) const override;

    ParamRefs & get_parameters(ParamRefs & params, ParamFilter * filter = nullptr) const override;
    ParamCRefs & get_parameters_const(ParamCRefs & params, ParamFilter * filter = nullptr) const override;

    void set_extra_param_factors(const Channel & channel, extra_param_factors & factors, size_t index
        ) const override;
    void set_grad_param_factors(const Channel & channel, grad_param_factors & factors, size_t index
        ) const override;
    
    std::string str() const override;

    Source(Components & components);
};

} // namespace fit
} // namespace gauss2d

#endif
