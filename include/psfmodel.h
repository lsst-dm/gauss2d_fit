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
    std::unique_ptr<gauss2d::Gaussians> get_gaussians(const Channel & channel = Channel::NONE()) const;

    ParamRefs & get_parameters(ParamRefs & params, ParamFilter * filter = nullptr) const override;
    ParamCRefs & get_parameters_const(ParamCRefs & params, ParamFilter * filter = nullptr) const override;

    std::string str() const override;

    PsfModel(Components & components);
    ~PsfModel();
};

} // namespace fit
} // namespace gauss2d

#endif