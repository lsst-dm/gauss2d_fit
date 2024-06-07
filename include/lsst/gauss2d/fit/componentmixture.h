#ifndef LSST_GAUSS2D_FIT_COMPONENTMIXTURE_H
#define LSST_GAUSS2D_FIT_COMPONENTMIXTURE_H

#include "component.h"
#include "parametricmodel.h"

namespace lsst::gauss2d::fit{

// Would consider making this a unique_ptr but pybind11 prevents this:
// e.g. https://github.com/pybind/pybind11/issues/1132
// and https://github.com/pybind/pybind11/issues/1161
// So the interface may as well support getting them

typedef std::vector<std::shared_ptr<Component>> Components;

/**
 * @brief A list of related Component instances.
 *
 * A ComponentMixture is any object that is composed of and can return a list
 * of its constituent components.
 */
class ComponentMixture : public ParametricModel {
public:
    virtual Components get_components() const = 0;
};

}  // namespace lsst::gauss2d::fit

#endif
