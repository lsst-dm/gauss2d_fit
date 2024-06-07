#include "lsst/gauss2d/type_name.h"

#include "lsst/gauss2d/fit/param_defs.h"

namespace lsst::gauss2d::fit {
double finite_difference_param(ParamBase& param, double delta) {
    const double value = param.get_value_transformed();
    double value_new = value + delta;
    try {
        param.set_value_transformed(value_new);
        // If the param is near an upper limit, this might fail: try -delta then
        if (param.get_value_transformed() != value_new) {
            delta = -delta;
            value_new = value + delta;
            param.set_value_transformed(value_new);
        }
    } catch (const std::runtime_error& err) {
        delta = -delta;
        value_new = value + delta;
        try {
            param.set_value_transformed(value_new);
        } catch (const std::runtime_error& err2) {
            throw std::runtime_error("Couldn't set param=" + param.str()
                                     + " to new value_transformed=" + std::to_string(value_new) + " due to "
                                     + err2.what() + "; are limits too restrictive?");
        }
    }
    return delta;
}
}  // namespace lsst::gauss2d::fit
