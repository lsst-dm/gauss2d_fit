#ifndef LSST_GAUSS2D_FIT_GAUSSIANMODELINTEGRAL_H
#define LSST_GAUSS2D_FIT_GAUSSIANMODELINTEGRAL_H

#include "channel.h"
#include "ellipticalcomponent.h"
#include "integralmodel.h"
#include "param_defs.h"
#include "param_filter.h"

namespace lsst::gauss2d::fit {

/**
 * A single-channel GaussianIntegral referencing a Parametric IntegralModel.
 */
class GaussianModelIntegral : public GaussianIntegral {
public:
    /**
     * Construct a GaussianModelIntegral instance for one Channel.
     *
     * @param channel The Channel for the integral.
     * @param integralmodel An IntegralModel valid for channel.
     */
    explicit GaussianModelIntegral(const Channel& channel,
                                   const std::shared_ptr<const IntegralModel> integralmodel);
    ~GaussianModelIntegral();

    double get_value() const override;
    void set_value(double value) override;

    std::string repr(bool name_keywords = false,
                     std::string_view namespace_separator = Object::CC_NAMESPACE_SEPARATOR) const override;
    std::string str() const override;

protected:
    const Channel& _channel;
    const std::shared_ptr<const IntegralModel> _integralmodel;
};

}  // namespace lsst::gauss2d::fit

#endif
