#ifndef LSST_GAUSS2D_FIT_MATH_H
#define LSST_GAUSS2D_FIT_MATH_H

#include <cmath>

namespace lsst::gauss2d::fit {

const double LOG_1 = log(1);
const double SQRT_2_PI = sqrt(2. * M_PI);
const double LOG_1_M_LOG_SQRT_2_PI = LOG_1 - log(SQRT_2_PI);

template <typename T>
T logit(T p) {
    return log(p / (1 + p));
}

template <typename T>
double logpdf_norm(T residual, T sigma) {
    return LOG_1_M_LOG_SQRT_2_PI - log(sigma) - residual * residual / 2;
}

template <template <typename...> class Container, class Value>
Value sum_iter(const Container<Value>& container) {
    Value sum = 0;
    for (const auto& value : container) sum += value;
    return sum;
}

}  // namespace lsst::gauss2d::fit

#endif  // LSST_GAUSS2D_FIT_MATH_H
