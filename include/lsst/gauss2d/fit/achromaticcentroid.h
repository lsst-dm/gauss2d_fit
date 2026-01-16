#ifndef LSST_GAUSS2D_FIT_ACHROMATICCENTROID_H
#define LSST_GAUSS2D_FIT_ACHROMATICCENTROID_H

#include <memory>
#include <ranges>

#include "centroidparameters.h"
#include "multichannelcentroid.h"

namespace lsst::gauss2d::fit {
/**
 * A Centroid with Parameters for x and y
 */
class AchromaticCentroid : public MultiChannelCentroid {
public:
    typedef std::map<std::reference_wrapper<const Channel>, std::shared_ptr<CentroidParameters>> Data;

    /**
     * Construct an AchromaticCentroid.
     *
     * @param centroid The CentroidParameters to return for every channel.
     */
    explicit AchromaticCentroid(std::shared_ptr<CentroidParameters> centroid = nullptr);

    std::shared_ptr<CentroidParameters> at(const Channel& channel) override;
    std::shared_ptr<const CentroidParameters> at(const Channel& channel) const override;

    std::shared_ptr<CentroidParameters> find(const Channel& channel) override;
    std::shared_ptr<const CentroidParameters> find(const Channel& channel) const override;

    std::vector<std::reference_wrapper<const Channel>> get_channels() const override;

    ParamRefs& get_parameters(ParamRefs& params, ParamFilter* filter = nullptr) const override;
    ParamCRefs& get_parameters_const(ParamCRefs& params, ParamFilter* filter = nullptr) const override;

    size_t n_channels() const override;

    std::string repr(bool name_keywords = false,
                     std::string_view namespace_separator = Object::CC_NAMESPACE_SEPARATOR) const override;
    std::string str() const override;

    ~AchromaticCentroid() override;

protected:
    std::shared_ptr<CentroidParameters> _centroid;
};
}  // namespace lsst::gauss2d::fit

#endif
