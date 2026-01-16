#include "lsst/gauss2d/type_name.h"

#include "lsst/gauss2d/fit/achromaticcentroid.h"

#include "lsst/gauss2d/fit/sersicparametricellipse.h"

namespace lsst::gauss2d::fit {

AchromaticCentroid::AchromaticCentroid(std::shared_ptr<CentroidParameters> centroid)
        : _centroid(centroid ? std::move(centroid) : std::make_shared<g2f::CentroidParameters>()) {}

std::shared_ptr<CentroidParameters> AchromaticCentroid::at(const Channel& channel) { return _centroid; }

std::shared_ptr<const CentroidParameters> AchromaticCentroid::at(const Channel& channel) const {
    return _centroid;
}

std::shared_ptr<CentroidParameters> AchromaticCentroid::find(const Channel& channel) { return _centroid; }

std::shared_ptr<const CentroidParameters> AchromaticCentroid::find(const Channel& channel) const {
    return _centroid;
}

std::vector<std::reference_wrapper<const Channel>> AchromaticCentroid::get_channels() const { return {}; }

ParamRefs& AchromaticCentroid::get_parameters(ParamRefs& params, ParamFilter* filter) const {
    return _centroid->get_parameters(params, filter);
}

ParamCRefs& AchromaticCentroid::get_parameters_const(ParamCRefs& params, ParamFilter* filter) const {
    return _centroid->get_parameters_const(params, filter);
}

size_t AchromaticCentroid::n_channels() const { return 0; }

std::string AchromaticCentroid::repr(bool name_keywords, std::string_view namespace_separator) const {
    std::string s = type_name_str<AchromaticCentroid>(false, namespace_separator) + "("
                    + (name_keywords ? "centroid=" : "") + _centroid->repr(name_keywords, namespace_separator)
                    + ")";
    return s;
}

std::string AchromaticCentroid::str() const {
    std::string s = type_name_str<AchromaticCentroid>(true) + "(centroid=" + _centroid->str() + ")";
    return s;
}

AchromaticCentroid::~AchromaticCentroid() = default;

}  // namespace lsst::gauss2d::fit