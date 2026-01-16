#include <ranges>
#include <stdexcept>

#include "lsst/gauss2d/fit/ellipticalcomponent.h"

#include "lsst/gauss2d/fit/achromaticcentroid.h"
#include "lsst/gauss2d/fit/centroidparameters.h"
#include "lsst/gauss2d/fit/component.h"
#include "lsst/gauss2d/fit/chromaticcentroid.h"
#include "lsst/gauss2d/fit/parametricellipse.h"
#include "lsst/gauss2d/fit/linearintegralmodel.h"
#include "lsst/gauss2d/fit/param_defs.h"
#include "lsst/gauss2d/fit/param_filter.h"

#include <unordered_set>

namespace lsst::gauss2d::fit {

EllipticalComponent::EllipticalComponent(std::shared_ptr<ParametricEllipse> ellipse,
                                         std::shared_ptr<MultiChannelCentroid> centroid,
                                         std::shared_ptr<IntegralModel> integralmodel)
        : _ellipse(std::move(ellipse)),
          _centroid(std::move(centroid)),
          _integralmodel(std::move(integralmodel)) {
    const bool has_integral = _integralmodel != nullptr;
    if (!has_integral) {
        _integralmodel = std::make_shared<LinearIntegralModel>(nullptr);
    }
    const auto channels_integral = _integralmodel->get_channels();
    std::unordered_set<std::string> channelset_integral;
    std::set<std::reference_wrapper<const Channel>> channelset_ordered_integral;
    const bool is_centroid_null = _centroid == nullptr;
    for (const auto& channel : channels_integral) {
        channelset_integral.insert(channel.get().name);
        if (is_centroid_null) {
            channelset_ordered_integral.insert(channel);
        }
    }
    if (is_centroid_null) {
        if (has_integral) {
            auto centroid_ptr = std::make_shared<CentroidParameters>();
            auto data = ChromaticCentroid::Data();
            for (const auto& channel : channelset_ordered_integral) {
                data.insert({channel, centroid_ptr});
            }
            _centroid = std::make_shared<ChromaticCentroid>(data);
        } else {
            _centroid = std::make_shared<AchromaticCentroid>();
        }
    } else {
        const auto channels_centroid = _centroid->get_channels();
        // If the channels are empty, then it's achromatic and there's no need
        // to check channel set compatibility
        if (!channels_centroid.empty()) {
            std::unordered_set<std::string> channelset_centroid;
            for (const auto& channel : channels_centroid) {
                channelset_centroid.insert(channel.get().name);
            }
            if (channelset_centroid != channelset_integral) {
                throw std::invalid_argument("centroid channels=" + to_string_iter(channelset_centroid)
                                            + " != integral channels=" + to_string_iter(channelset_integral));
            }
        }
    }
    // Superclasses are responsible for default-constructing ellipse
    if (_ellipse == nullptr) throw std::invalid_argument("ellipse must not be null");
}

const MultiChannelCentroid& EllipticalComponent::get_centroid() const { return *_centroid; }
std::vector<std::reference_wrapper<const Channel>> EllipticalComponent::get_channels() const {
    return _centroid->get_channels();
}
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
