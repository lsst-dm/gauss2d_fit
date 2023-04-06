#include "ellipticalcomponent.h"

#include "centroidparameters.h"
#include "component.h"
#include "parametricellipse.h"
#include "linearintegralmodel.h"
#include "param_defs.h"
#include "param_filter.h"
#include <stdexcept>

namespace gauss2d::fit {

const CentroidParameters& EllipticalComponent::get_centroid() const { return *_centroid; }
const ParametricEllipse& EllipticalComponent::get_ellipse() const { return *_ellipse; }
const IntegralModel& EllipticalComponent::get_integralmodel() const { return *_integralmodel; }

ParamRefs& EllipticalComponent::get_parameters(ParamRefs& params, ParamFilter* filter) const {
    _centroid->get_parameters(params, filter);
    _ellipse->get_parameters(params, filter);
    _integralmodel->get_parameters(params, filter);
    return params;
}
ParamCRefs& EllipticalComponent::get_parameters_const(ParamCRefs& params, ParamFilter* filter) const {
    _centroid->get_parameters_const(params, filter);
    _ellipse->get_parameters_const(params, filter);
    _integralmodel->get_parameters_const(params, filter);
    return params;
}

std::string EllipticalComponent::repr(bool name_keywords) const {
    return (name_keywords ? "centroid=" : "") + _centroid->repr(name_keywords) + ", "
           + (name_keywords ? "ellipse=" : "") + _ellipse->repr(name_keywords) + ", "
           + (name_keywords ? "integralmodel=" : "") + _integralmodel->repr(name_keywords);
}

std::string EllipticalComponent::str() const {
    return "centroid=" + _centroid->str() + ", ellipse=" + _ellipse->str()
           + ", integralmodel=" + _integralmodel->str();
}

EllipticalComponent::EllipticalComponent(std::shared_ptr<ParametricEllipse> ellipse,
                                         std::shared_ptr<CentroidParameters> centroid,
                                         std::shared_ptr<IntegralModel> integralmodel)
        : _ellipse(std::move(ellipse)),
          _centroid(centroid ? std::move(centroid) : std::make_shared<CentroidParameters>()),
          _integralmodel(integralmodel ? std::move(integralmodel)
                                       : std::make_shared<LinearIntegralModel>(nullptr)) {
    if (_ellipse == nullptr) throw std::invalid_argument("ellipse must not be null");
}

}  // namespace gauss2d::fit
