#ifndef GAUSS2D_FIT_PARAMETRICMODEL_H
#define GAUSS2D_FIT_PARAMETRICMODEL_H

#include "gauss2d/gaussian.h"

#include "channel.h"
#include "parametric.h"

namespace gauss2d
{
namespace fit
{

class ParametricModel : public Parametric
{
public:
    virtual std::unique_ptr<const gauss2d::Gaussians> get_gaussians(const Channel & channel) const = 0;
};

} // namespace fit
} // namespace gauss2d

#endif