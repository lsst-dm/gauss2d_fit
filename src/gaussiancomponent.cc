#include "gaussiancomponent.h"

#include "gauss2d/gaussian.h"
#include "integralmodel.h"
#include "gaussianparametricellipse.h"
#include "gaussianmodelintegral.h"
#include "param_defs.h"
#include "param_filter.h"
#include <stdexcept>

namespace gauss2d
{
namespace fit
{

std::unique_ptr<const gauss2d::Gaussians> GaussianComponent::get_gaussians(const Channel & channel) const {
    gauss2d::Gaussians::Data gaussians = {
        std::make_shared<Gaussian>(
            std::make_shared<Centroid>(this->_centroid),
            std::make_shared<Ellipse>(this->_ellipsedata),
            std::make_shared<GaussianModelIntegral>(channel, this->_integralmodel)
        )
    };
    return std::make_unique<const gauss2d::Gaussians>(gaussians);
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
    std::shared_ptr<GaussianParametricEllipse> ellipse,
    std::shared_ptr<CentroidParameters> centroid,
    std::shared_ptr<IntegralModel> integralmodel
) : GaussianParametricEllipseHolder(std::move(ellipse)), EllipticalComponent(_ellipsedata, centroid, integralmodel)
{
}

} // namespace fit
} // namespace gauss2d
