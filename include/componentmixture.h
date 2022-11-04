#ifndef GAUSS2D_FIT_COMPONENTMIXTURE_H
#define GAUSS2D_FIT_COMPONENTMIXTURE_H

#include "component.h"
#include "parametricmodel.h"

namespace gauss2d
{
namespace fit
{

// Would consider making this a unique_ptr but pybind11 prevents this:
// e.g. https://github.com/pybind/pybind11/issues/1132
// and https://github.com/pybind/pybind11/issues/1161
// So the interface may as well support getting them

typedef std::vector<std::shared_ptr<Component>> Components;

class ComponentMixture : public ParametricModel
{
public:
    virtual Components get_components() const = 0;
};

} // namespace fit
} // namespace gauss2d

#endif
