#include "param_defs.h"

namespace gauss2d::fit {
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
        param.set_value_transformed(value_new);
        if (param.get_value_transformed() != value_new) {
            throw std::runtime_error("Couldn't set param=" + param.str()
                                     + " to new value_transformed=" + std::to_string(value) + " +/- "
                                     + std::to_string(delta) + "; are limits too restrictive?");
        }
    }
    return delta;
}
}  // namespace gauss2d::fit
