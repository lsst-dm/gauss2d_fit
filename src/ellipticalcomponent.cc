#include "ellipticalcomponent.h"

#include "centroidparameters.h"
#include "component.h"
#include "ellipseparameters.h"
#include "linearintegralmodel.h"
#include "param_defs.h"
#include "param_filter.h"

namespace gauss2d
{
namespace fit
{

const CentroidParameters & EllipticalComponent::get_centroid() const { return *_centroid; }
const EllipseParameters & EllipticalComponent::get_ellipse() const { return *_ellipse; }
const IntegralModel & EllipticalComponent::get_integralmodel() const { return *_integralmodel; }

ParamRefs & EllipticalComponent::get_parameters(ParamRefs & params, ParamFilter * filter) const {
    _centroid->get_parameters(params, filter);
    _ellipse->get_parameters(params, filter);
    _integralmodel->get_parameters(params, filter);
    return params;
}
ParamCRefs & EllipticalComponent::get_parameters_const(ParamCRefs & params, ParamFilter * filter) const {
    _centroid->get_parameters_const(params, filter);
    _ellipse->get_parameters_const(params, filter);
    _integralmodel->get_parameters_const(params, filter);
    return params;

}

std::string EllipticalComponent::str() const {
    return "centroid=" + _centroid->str() + ", ellipse=" + _ellipse->str()
        + ", integralmodel=" + _integralmodel->str();
}

EllipticalComponent::EllipticalComponent(
    std::shared_ptr<CentroidParameters> centroid,
    std::shared_ptr<EllipseParameters> ellipse,
    std::shared_ptr<IntegralModel> integralmodel
) : _centroid(centroid ? std::move(centroid) : std::make_shared<CentroidParameters>()),
    _ellipse(ellipse ? std::move(ellipse) : std::make_shared<EllipseParameters>()),
    _integralmodel(integralmodel ? std::move(integralmodel) : 
        std::make_shared<LinearIntegralModel>(nullptr))
{}

} // namespace fit
} // namespace gauss2d
