#ifndef GAUSS2D_FIT_CHROMATIC_H
#define GAUSS2D_FIT_CHROMATIC_H

#include "channel.h"

namespace gauss2d::fit {

class Chromatic {
public:
    /// @brief Get the set of channels this instance is applicable for
    /// @note Implementers must return a set (all unique items). This cannot be enforced (yet).
    virtual std::vector<std::reference_wrapper<const Channel>> get_channels() const = 0;

    virtual ~Chromatic() = default;
};

}  // namespace gauss2d::fit

#endif
