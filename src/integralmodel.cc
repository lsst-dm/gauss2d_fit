#ifndef GAUSS2D_FIT_INTEGRALMODEL_H
#define GAUSS2D_FIT_INTEGRALMODEL_H

#include <memory>

#include "channel.h"
#include "parametric.h"

namespace gauss2d::fit {

class IntegralModel : public Parametric {
public:
    virtual double get_integral(const Channel& channel) const = 0;
};

}  // namespace gauss2d::fit

#endif