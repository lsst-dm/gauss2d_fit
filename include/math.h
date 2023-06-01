#ifndef GAUSS2DFIT_MATH_H
#define GAUSS2DFIT_MATH_H

#include <cmath>

namespace gauss2d::fit {

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

}  // namespace gauss2d::fit

#endif  // GAUSS2DFIT_MATH_H
