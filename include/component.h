#ifndef GAUSS2D_FIT_COMPONENT_H
#define GAUSS2D_FIT_COMPONENT_H

#include "channel.h"
#include "integralmodel.h"
#include "parametricmodel.h"

namespace gauss2d::fit
{

class Component : public ParametricModel
{
public:
    virtual const IntegralModel & get_integralmodel() const = 0;
};

} // namespace gauss2d::fit

#endif