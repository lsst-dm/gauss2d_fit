#ifndef LSST_GAUSS2D_FIT_PARAMETRIC_H
#define LSST_GAUSS2D_FIT_PARAMETRIC_H

#include "lsst/gauss2d/object.h"

#include "param_defs.h"
#include "param_filter.h"

namespace lsst::gauss2d::fit{
/**
 * A parametric object that can return and filter its Parameter instances.
 */
class Parametric : public Object {
public:
    /**
     * Add Parameter refs matching the filter to a vector, in order.
     *
     * @param params The vector to add to.
     * @param filter The filter to apply to this Object's parameters.
     * @return A ref to params (for method chaining)
     */
    virtual ParamRefs& get_parameters(ParamRefs& params, ParamFilter* filter = nullptr) const = 0;
    /// Same as get_parameters(), but for const refs.
    virtual ParamCRefs& get_parameters_const(ParamCRefs& params, ParamFilter* filter = nullptr) const = 0;

    /// Same as get_parameters(), but returning a new vector.
    ParamRefs get_parameters_new(ParamFilter* filter = nullptr) const {
        ParamRefs params{};
        get_parameters(params, filter);
        return params;
    }
    /// Same as get_parameters_const(), but returning a new vector.
    ParamCRefs get_parameters_const_new(ParamFilter* filter = nullptr) const {
        ParamCRefs params{};
        get_parameters_const(params, filter);
        return params;
    }

    virtual ~Parametric() = default;
};
}  // namespace lsst::gauss2d::fit

#endif
