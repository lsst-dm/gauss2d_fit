#include "lsst/gauss2d/type_name.h"

#include "lsst/gauss2d/fit/chromaticcentroid.h"

#include "lsst/gauss2d/fit/sersicparametricellipse.h"

namespace lsst::gauss2d::fit {

ChromaticCentroid::ChromaticCentroid(std::optional<Data> data) {
    if (data == std::nullopt) {
        _data[Channel::NONE()] = std::make_shared<g2f::CentroidParameters>();
    } else {
        for (auto [channel, centroid] : *data) {
            _data[channel] = centroid;
        }
    }
}

std::shared_ptr<CentroidParameters> ChromaticCentroid::at(const Channel& channel) {
    return _data.at(channel);
}

std::shared_ptr<const CentroidParameters> ChromaticCentroid::at(const Channel& channel) const {
    return _data.at(channel);
}

std::shared_ptr<CentroidParameters> ChromaticCentroid::find(const Channel& channel) {
    auto found = _data.find(channel);
    return found == _data.end() ? nullptr : found->second;
}

std::shared_ptr<const CentroidParameters> ChromaticCentroid::find(const Channel& channel) const {
    auto found = _data.find(channel);
    return found == _data.end() ? nullptr : found->second;
}

std::vector<std::reference_wrapper<const Channel>> ChromaticCentroid::get_channels() const {
    std::vector<std::reference_wrapper<const Channel>> rval = {};
    for (const auto& channel : _data | std::views::keys) {
        rval.emplace_back(channel);
    }
    return rval;
}

template <typename P, bool get_const>
P& _get_parameters(P& params, ParamFilter* filter, const ChromaticCentroid::Data& data) {
    const auto channel_filter
            = ((filter != nullptr) && (filter->channel != std::nullopt)) ? filter->channel : std::nullopt;
    const bool has_channel = channel_filter != std::nullopt;
    for (const auto& [channel, centroid] : data) {
        if (!has_channel || (channel_filter == channel)) {
            if constexpr (get_const) {
                centroid->get_parameters_const(params, filter);
            } else {
                centroid->get_parameters(params, filter);
            }
        }
    }
    return params;
}

ParamRefs& ChromaticCentroid::get_parameters(ParamRefs& params, ParamFilter* filter) const {
    return _get_parameters<ParamRefs, false>(params, filter, _data);
}

ParamCRefs& ChromaticCentroid::get_parameters_const(ParamCRefs& params, ParamFilter* filter) const {
    return _get_parameters<ParamCRefs, true>(params, filter, _data);
}

size_t ChromaticCentroid::n_channels() const { return _data.size(); }

std::string ChromaticCentroid::repr(bool name_keywords, std::string_view namespace_separator) const {
    std::string s = type_name_str<ChromaticCentroid>(false, namespace_separator) + "("
                    + (name_keywords ? "data=" : "") + "{";
    for (auto [channel, centroid] : _data) {
        s += channel.get().repr(name_keywords, namespace_separator) + ": "
             + centroid->repr(name_keywords, namespace_separator) + ",";
    }
    return s + "})";
}

std::string ChromaticCentroid::str() const {
    std::string s = type_name_str<ChromaticCentroid>(true) + "(data={";
    for (const auto& [channel, centroid] : _data) {
        s += channel.get().str() + ": " + centroid->str() + ",";
    }
    return s + "})";
}

ChromaticCentroid::~ChromaticCentroid() = default;

}  // namespace lsst::gauss2d::fit