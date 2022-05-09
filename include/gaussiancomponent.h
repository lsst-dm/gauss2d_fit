#ifndef GAUSS2D_FIT_GAUSSIANCOMPONENT_H
#define GAUSS2D_FIT_GAUSSIANCOMPONENT_H

#include "component.h"
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
private:
    const Channel & _channel;
    const std::shared_ptr<const IntegralModel> _integralmodel;

public:
    double get_value() const;
    void set_value(double value);

    std::string str() const override;

    GaussianModelIntegral(
        const Channel & channel, const std::shared_ptr<const IntegralModel> integralmodel
    );
    ~GaussianModelIntegral();
};

class GaussianComponent : public EllipticalComponent {
public:
    std::unique_ptr<gauss2d::Gaussians> get_gaussians(const Channel & channel) const;

    ParamRefs & get_parameters(ParamRefs & params, ParamFilter * filter = nullptr) const override;
    ParamCRefs & get_parameters_const(ParamCRefs & params, ParamFilter * filter = nullptr) const override;

    std::string str() const override;
    
    GaussianComponent(
        std::shared_ptr<CentroidParameters> centroid = nullptr,
        std::shared_ptr<EllipseParameters> ellipse = nullptr,
        std::shared_ptr<IntegralModel> integralmodel = nullptr
    );
};
} // namespace fit
} // namespace gauss2d

#endif