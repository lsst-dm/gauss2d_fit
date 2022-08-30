#ifndef GAUSS2D_FIT_PSFMODEL_H
#define GAUSS2D_FIT_PSFMODEL_H

#include <memory>

#include "component.h"
#include "param_filter.h"
#include "parametricmodel.h"

namespace gauss2d
{
namespace fit
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
class PsfModel : public ParametricModel
{
public:
    // Would like this to be unique_ptr but can't due to various pybind issues
    // e.g. https://github.com/pybind/pybind11/issues/1132
    // and https://github.com/pybind/pybind11/issues/1161
    typedef std::vector<std::shared_ptr<Component>> Components;
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
    
    std::unique_ptr<const gauss2d::Gaussians> get_gaussians(const Channel & channel = Channel::NONE()) const override;

    ParamRefs & get_parameters(ParamRefs & params, ParamFilter * filter = nullptr) const override;
    ParamCRefs & get_parameters_const(ParamCRefs & params, ParamFilter * filter = nullptr) const override;

    std::string str() const override;

    PsfModel(Components & components);
    ~PsfModel();
};

} // namespace fit
} // namespace gauss2d

#endif
