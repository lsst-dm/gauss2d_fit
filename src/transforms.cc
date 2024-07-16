#include <limits>

#include "lsst/gauss2d/type_name.h"

#include "lsst/gauss2d/fit/transforms.h"

namespace lsst::gauss2d::fit {
static const double INF = std::numeric_limits<double>::infinity();

double LogitLimitedTransform::derivative(double x) const {
    double y = (x - _limits->get_min()) / _range;
    if (y == 1)
        return INF;
    else if (y == 0)
        return -INF;
    return (1 / y + 1 / (1 - y)) * _factor / _range;
}

double LogitLimitedTransform::forward(double x) const {
    double min = _limits->get_min();
    double max = _limits->get_max();
    if (x == min)
        return -std::numeric_limits<double>::infinity();
    else if (x == max)
        return std::numeric_limits<double>::infinity();
    double y = (x - min) / _range;
    if (!(y < 1) || !(y > 0)) return nan("");
    return log(y / (1 - y)) * _factor;
}

double LogitLimitedTransform::reverse(double x) const {
    double y = -x * _factor;
    // math.log(np.finfo(np.float64) = 709.782712893384
    // things will go badly well before then
    if (y > 709.7827) return _limits->get_min();
    y = 1 + exp(y);
    return _range / y + _limits->get_min();
}
}  // namespace lsst::gauss2d::fit
