#ifndef GAUSS2D_FIT_PARAM_DEFS_H
#define GAUSS2D_FIT_PARAM_DEFS_H

#include <vector>

#include "parameters/parameter.h"

namespace gauss2d::fit {
using ParamBase = parameters::ParameterBase<double>;
using ParamBaseCRef = std::reference_wrapper<const ParamBase>;
using ParamCRefs = std::vector<ParamBaseCRef>;
using ParamBaseRef = std::reference_wrapper<ParamBase>;
using ParamRefs = std::vector<ParamBaseRef>;
}  // namespace gauss2d::fit
#endif  // GAUSS2D_FIT_PARAM_DEFS_H
