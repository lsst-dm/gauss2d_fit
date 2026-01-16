#ifndef LSST_GAUSS2D_FIT_MULTICHANNELCENTROID_H
#define LSST_GAUSS2D_FIT_MULTICHANNELCENTROID_H

#include <memory>

#include "centroidparameters.h"
#include "chromatic.h"
#include "parametric.h"

namespace lsst::gauss2d::fit {
/**
 * An interface made to return CentroidParameters by Channel.
 */
class MultiChannelCentroid : public Chromatic, public Parametric {
public:
    /// Get the CentroidParameters for the given Channel
    virtual std::shared_ptr<CentroidParameters> at(const Channel& channel) = 0;
    /// Get the (const) CentroidParameters for the given Channel
    virtual std::shared_ptr<const CentroidParameters> at(const Channel& channel) const = 0;

    /// Find the CentroidParameters for the given Channel
    virtual std::shared_ptr<CentroidParameters> find(const Channel& channel) = 0;
    /// Find the (const) CentroidParameters for the given Channel
    virtual std::shared_ptr<const CentroidParameters> find(const Channel& channel) const = 0;

    /// Return the number of Channels supported, or 0 if any Channel is valid.
    virtual size_t n_channels() const = 0;

    virtual ~MultiChannelCentroid() {};
};
}  // namespace lsst::gauss2d::fit

#endif
