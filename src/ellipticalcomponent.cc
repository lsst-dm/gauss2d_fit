#include <stdexcept>

#include "lsst/gauss2d/fit/ellipticalcomponent.h"
#include "lsst/gauss2d/fit/centroidparameters.h"
#include "lsst/gauss2d/fit/component.h"
#include "lsst/gauss2d/fit/parametricellipse.h"
#include "lsst/gauss2d/fit/linearintegralmodel.h"
#include "lsst/gauss2d/fit/param_defs.h"
#include "lsst/gauss2d/fit/param_filter.h"

namespace lsst::gauss2d::fit {

EllipticalComponent::EllipticalComponent(std::shared_ptr<ParametricEllipse> ellipse,
                                         std::shared_ptr<CentroidParameters> centroid,
                                         std::shared_ptr<IntegralModel> integralmodel)
        : _ellipse(std::move(ellipse)),
          _centroid(centroid ? std::move(centroid) : std::make_shared<CentroidParameters>()),
          _integralmodel(integralmodel ? std::move(integralmodel)
                                       : std::make_shared<LinearIntegralModel>(nullptr)) {
    if (_ellipse == nullptr) throw std::invalid_argument("ellipse must not be null");
}

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

std::string EllipticalComponent::repr(bool name_keywords, std::string_view namespace_separator) const {
    return (name_keywords ? "ellipse=" : "") + _ellipse->repr(name_keywords, namespace_separator) + ", "
           + (name_keywords ? "centroid=" : "") + _centroid->repr(name_keywords, namespace_separator) + ", "
           + (name_keywords ? "integralmodel=" : "")
           + _integralmodel->repr(name_keywords, namespace_separator);
}

std::string EllipticalComponent::str() const {
    return "ellipse=" + _ellipse->str() + ", centroid=" + _centroid->str()
           + ", integralmodel=" + _integralmodel->str();
}

}  // namespace lsst::gauss2d::fit
