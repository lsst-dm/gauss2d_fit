#ifndef LSST_GAUSS2D_FIT_PARAM_DEFS_H
#define LSST_GAUSS2D_FIT_PARAM_DEFS_H

#include <vector>

#include "lsst/modelfit/parameters.h"

namespace lsst::gauss2d::fit{
using ParamBase = lsst::modelfit::parameters::ParameterBase<double>;
using ParamBaseCRef = std::reference_wrapper<const ParamBase>;
using ParamCRefs = std::vector<ParamBaseCRef>;
using ParamBaseRef = std::reference_wrapper<ParamBase>;
using ParamRefs = std::vector<ParamBaseRef>;

double finite_difference_param(ParamBase& param, double delta);
/**
 * Increment a parameter by a small difference delta.
 *
 * @param param The parameter to increment.
 * @param delta The target delta value.
 * @return The actual delta value.
 *
 * @note The delta value will switch sign if the target fails.
 */

}  // namespace lsst::gauss2d::fit
#endif  // LSST_GAUSS2D_FIT_PARAM_DEFS_H
