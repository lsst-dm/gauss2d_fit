#ifndef GAUSS2D_FIT_CENTROIDPARAMETERS_H
#define GAUSS2D_FIT_CENTROIDPARAMETERS_H

#include <memory>

#include "gauss2d/centroid.h"
#include "parameters.h"

namespace gauss2d
{
namespace fit
{

class CentroidParameters : public gauss2d::CentroidData
{
private:
    std::shared_ptr<CentroidXParameter> _x;
    std::shared_ptr<CentroidYParameter> _y;

public:
    double get_x() const { return _x->get_value(); }
    double get_y() const { return _y->get_value(); }
    std::array<double, 2> get_xy() const { return {this->get_x(), this->get_y()}; }

    CentroidXParameter & get_x_param() const { return *_x; }
    CentroidYParameter & get_y_param() const { return *_y; }

    std::shared_ptr<CentroidXParameter> get_x_param_ptr() { return _x; }
    std::shared_ptr<CentroidYParameter> get_y_param_ptr() { return _y; }

    void set_x(double x) { _x->set_value(x); }
    void set_y(double y) { _y->set_value(y); }
    void set_xy(const std::array<double, 2> & xy) { this->set_x(xy[0]); this->set_y(xy[1]); }

    std::string str() const {
        return "CentroidParameters(x=" + _x->str() + ", y=" + _y->str() + ")";
    }

    CentroidParameters(
        std::shared_ptr<CentroidXParameter> x = nullptr,
        std::shared_ptr<CentroidYParameter> y = nullptr
    ) :
        CentroidData(),
        _x(x == nullptr ? std::make_shared<CentroidXParameter>() : std::move(x)),
        _y(y == nullptr ? std::make_shared<CentroidYParameter>() : std::move(y))
    {}

    CentroidParameters(double x, double y) : CentroidData(), 
        _x(std::make_shared<CentroidXParameter>(x)),
        _y(std::make_shared<CentroidYParameter>(y))
    {}
};
} // namespace fit
} // namespace gauss2d

#endif