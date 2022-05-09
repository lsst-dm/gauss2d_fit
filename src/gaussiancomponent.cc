#include "gaussiancomponent.h"

#include <iostream>
#include <stdexcept>

#include "component.h"
#include "ellipticalcomponent.h"
#include "gauss2d/gaussian.h"
#include "integralmodel.h"
#include "linearintegralmodel.h"
#include "param_defs.h"
#include "param_filter.h"

namespace gauss2d
{
namespace fit
{

double GaussianModelIntegral::get_value() const {
    return _integralmodel->get_integral(_channel);
}
void GaussianModelIntegral::set_value(double value) {
    throw std::runtime_error("Can't set_value on GaussianModelIntegral");
}

std::string GaussianModelIntegral::str() const {
    return "GaussianModelIntegral(channel=" + _channel.str()
        + ", integralmodel=" + _integralmodel->str() + ")";
}

GaussianModelIntegral::GaussianModelIntegral(
    const Channel & channel, const std::shared_ptr<const IntegralModel> integralmodel
) : _channel(channel), _integralmodel(std::move(integralmodel))
{
    if(_integralmodel == nullptr) throw std::invalid_argument("GaussianModelIntegral integralmodel can't be null");
}
GaussianModelIntegral::~GaussianModelIntegral() {};

std::unique_ptr<gauss2d::Gaussians> GaussianComponent::get_gaussians(const Channel & channel) const {
    gauss2d::Gaussians::Data gaussians = {
        std::make_shared<Gaussian>(
            std::make_shared<Centroid>(this->_centroid),
            std::make_shared<Ellipse>(this->_ellipse),
            std::make_shared<GaussianModelIntegral>(channel, this->_integralmodel)
        )
    };
    return std::make_unique<gauss2d::Gaussians>(gaussians);
}

ParamRefs & GaussianComponent::get_parameters(ParamRefs & params, ParamFilter * filter) const {
    EllipticalComponent::get_parameters(params, filter);
    return params;
}

ParamCRefs & GaussianComponent::get_parameters_const(ParamCRefs & params, ParamFilter * filter) const {
    EllipticalComponent::get_parameters_const(params, filter);
    return params;
}

std::string GaussianComponent::str() const {
    return "GaussianComponent(" + EllipticalComponent::str() + ", photo=" + _integralmodel->str() + ")";
}

GaussianComponent::GaussianComponent(
    std::shared_ptr<CentroidParameters> centroid,
    std::shared_ptr<EllipseParameters> ellipse,
    std::shared_ptr<IntegralModel> integralmodel
) : EllipticalComponent(centroid, ellipse, integralmodel)
{
}

} // namespace fit
} // namespace gauss2d
