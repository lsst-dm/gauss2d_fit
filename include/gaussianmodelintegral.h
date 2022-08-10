#ifndef GAUSS2D_FIT_GAUSSIANMODELINTEGRAL_H
#define GAUSS2D_FIT_GAUSSIANMODELINTEGRAL_H

#include "channel.h"
#include "ellipticalcomponent.h"
#include "integralmodel.h"
#include "param_defs.h"
#include "param_filter.h"

namespace gauss2d
{
namespace fit
{

class GaussianModelIntegral : public GaussianIntegral
{
protected:
    const Channel & _channel;
    const std::shared_ptr<const IntegralModel> _integralmodel;

public:
    virtual double get_value() const override;
    virtual void set_value(double value) override;

    virtual std::string str() const override;

    GaussianModelIntegral(
        const Channel & channel, const std::shared_ptr<const IntegralModel> integralmodel
    );
    ~GaussianModelIntegral();
};

} // namespace fit
} // namespace gauss2d

#endif