#ifndef GAUSS2D_FIT_SOURCE_H
#define GAUSS2D_FIT_SOURCE_H

#include <memory>

#include "component.h"
#include "parametricmodel.h"
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
class Source : public ParametricModel
{
public:
    // Would like this to be unique_ptr but can't due to various pybind issues
    // e.g. https://github.com/pybind/pybind11/issues/1132
    // and https://github.com/pybind/pybind11/issues/1161

    typedef std::vector<std::shared_ptr<Component>> Components;

private:
    Components _components = {};

public:
    std::unique_ptr<gauss2d::Gaussians> get_gaussians(const Channel & channel) const override;

    ParamRefs & get_parameters(ParamRefs & params, ParamFilter * filter = nullptr) const override;
    ParamCRefs & get_parameters_const(ParamCRefs & params, ParamFilter * filter = nullptr) const override;

    std::string str() const override;

    Source(Components & components);
};

} // namespace fit
} // namespace gauss2d

#endif