#ifndef GAUSS2D_FIT_GAUSSIANMODELINTEGRAL_H
#define GAUSS2D_FIT_GAUSSIANMODELINTEGRAL_H

#include "channel.h"
#include "ellipticalcomponent.h"
#include "integralmodel.h"
#include "param_defs.h"
#include "param_filter.h"

namespace gauss2d::fit
{

/**
 * A single-channel GaussianIntegral referencing a Parametric IntegralModel.
 */
class GaussianModelIntegral : public GaussianIntegral
{
protected:
    const Channel & _channel;
    const std::shared_ptr<const IntegralModel> _integralmodel;

public:
    double get_value() const override;
    void set_value(double value) override;

    std::string repr(bool name_keywords = false) const override;
    std::string str() const override;

    /**
     * Construct a GaussianModelIntegral instance for one Channel.
     *
     * @param channel The Channel for the integral.
     * @param integralmodel An IntegralModel valid for channel.
     */
    GaussianModelIntegral(
        const Channel & channel, const std::shared_ptr<const IntegralModel> integralmodel
    );
    ~GaussianModelIntegral();
};

} // namespace gauss2d::fit

#endif
