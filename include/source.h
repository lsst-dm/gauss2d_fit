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

class Source : public ParametricModel
{
public:
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