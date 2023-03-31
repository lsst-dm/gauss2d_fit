#ifndef GAUSS2D_FIT_PSFMODEL_H
#define GAUSS2D_FIT_PSFMODEL_H

#include <memory>

#include "component.h"
#include "componentmixture.h"
#include "param_filter.h"

namespace gauss2d::fit
{

/*
    Like a Source, a PsfModel is a collection of Components. However, a
    PsfModel is meant to represent the smoothing kernel for a single
    Observation (whether from the optical system, environmental conditions,
    or any other source of blurring). As such, it should have IntegralModels
    that sum to unity. This is most easily enforced with the use of
    FractionalIntegralModels.

    PsfModels are also generally required to not have a specific Channel.
    Logically, it should have the same Channel as the Observation it applies
    to, but generally, it cannot be defined to apply to multiple
    Observations, so not allowing non-NONE Channels is a compromise to reflect
    this fact.
*/
class PsfModel : public ComponentMixture
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
    void add_grad_param_factors(const Channel & channel, grad_param_factors & factor) const override;
    
    Components get_components() const override;
    std::unique_ptr<const gauss2d::Gaussians> get_gaussians(
        const Channel & channel = Channel::NONE()) const override;
    size_t get_n_gaussians(const Channel & channel = Channel::NONE()) const override;

    ParamRefs & get_parameters(ParamRefs & params, ParamFilter * filter = nullptr) const override;
    ParamCRefs & get_parameters_const(ParamCRefs & params, ParamFilter * filter = nullptr) const override;

    void set_extra_param_factors(const Channel & channel, extra_param_factors & factors, size_t index) const override;
    void set_grad_param_factors(const Channel & channel, grad_param_factors & factor, size_t index) const override;

    std::string repr(bool name_keywords = false) const override;
    std::string str() const override;

    explicit PsfModel(Components & components);
    ~PsfModel();
};

} // namespace gauss2d::fit

#endif
