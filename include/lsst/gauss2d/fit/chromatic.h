#ifndef LSST_GAUSS2D_FIT_CHROMATIC_H
#define LSST_GAUSS2D_FIT_CHROMATIC_H

#include "channel.h"

namespace lsst::gauss2d::fit{

class Chromatic {
public:
    /// @brief Get the set of channels this instance is applicable for
    /// @note Implementers must return a set (all unique items). This cannot be enforced (yet).
    virtual std::vector<std::reference_wrapper<const Channel>> get_channels() const = 0;

    virtual ~Chromatic() = default;
};

}  // namespace lsst::gauss2d::fit

#endif
