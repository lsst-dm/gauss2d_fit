#ifndef GAUSS2D_FIT_PARAMETRICELLIPSE_H
#define GAUSS2D_FIT_PARAMETRICELLIPSE_H

#include <array>

#include "parameters.h"
#include "parametric.h"
#include "param_defs.h"
#include "param_filter.h"

namespace gauss2d
{
namespace fit
{

class QuasiParametricEllipse : public Parametric
{
public:
    virtual double get_rho() const = 0;
    virtual double get_size_x() const = 0;
    virtual double get_size_y() const = 0;
    virtual std::array<double, 3> get_xyr() const { return {get_size_x(), get_size_y(), get_rho()}; };
};

class ParametricEllipse : public QuasiParametricEllipse
{
public:
    virtual RhoParameter & get_rho_param() const = 0;
    virtual SizeXParameter & get_size_x_param() const = 0;
    virtual SizeYParameter & get_size_y_param() const = 0;

    virtual std::shared_ptr<RhoParameter> get_rho_param_ptr() = 0;
    virtual std::shared_ptr<SizeXParameter> get_size_x_param_ptr() = 0;
    virtual std::shared_ptr<SizeYParameter> get_size_y_param_ptr() = 0;

    // No set because it would be ambiguous with gauss2d::set
    // TODO: Consider disambiguating
    virtual void set_rho(double rho) = 0;
    virtual void set_size_x(double size_x) = 0;
    virtual void set_size_y(double size_y) = 0;
    virtual void set_xyr(const std::array<double, 3> & xyr) {
        set_size_x(xyr[0]);
        set_size_y(xyr[1]);
        set_rho(xyr[2]);
    };
};

} // namespace fit
} // namespace gauss2d

#endif
