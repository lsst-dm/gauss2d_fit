#ifndef GAUSS2D_FIT_ELLIPSEPARAMETERS_H
#include "ellipseparameters.h"
#endif

namespace gauss2d
{
namespace fit
{
    std::array<double, 3> EllipseParameters::get_xyr() const {
        return {this->get_sigma_x(), this->get_sigma_y(), this->get_rho()};
    }

    void EllipseParameters::set(double sigma_x, double sigma_y, double rho) {
        set_sigma_x(sigma_x);
        set_sigma_y(sigma_y);
        set_rho(rho);
    }
} // namespace fit
} // namespace gauss2d
