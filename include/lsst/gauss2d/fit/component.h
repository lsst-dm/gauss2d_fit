#ifndef LSST_GAUSS2D_FIT_COMPONENT_H
#define LSST_GAUSS2D_FIT_COMPONENT_H

#include "channel.h"
#include "integralmodel.h"
#include "parametricmodel.h"

namespace lsst::gauss2d::fit{
/**
 * @brief An atomic constituent of a source.
 *
 * A Component is the smallest (atomic) constituent of a source, representing some 2D intensity distribution.
 * It can be composed of multiple Gaussians but must control all of their values with Parameters.
 */
class Component : public ParametricModel {
public:
    virtual const IntegralModel& get_integralmodel() const = 0;
};

}  // namespace lsst::gauss2d::fit

#endif