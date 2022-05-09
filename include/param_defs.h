#ifndef GAUSS2D_FIT_PARAM_DEFS_H
#define GAUSS2D_FIT_PARAM_DEFS_H

#include <vector>

#include "parameters/parameter.h"

namespace gauss2d
{
namespace fit
{
using ParamBase = parameters::ParameterBase<double>;
using ParamBaseCRef = std::reference_wrapper<const ParamBase>;
using ParamCRefs = std::vector<ParamBaseCRef>;
using ParamBaseRef = std::reference_wrapper<ParamBase>;
using ParamRefs = std::vector<ParamBaseRef>;
}
}
#endif //GAUSS2DFIT_PARAM_DEFS_H
