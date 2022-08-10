#include "centroidparameters.h"

#include <memory>

#include "parameters.h"

namespace gauss2d
{
namespace fit
{
template<typename t>
t & _get_parameters(
    t & params, ParamFilter * filter,
    CentroidXParameter & x, CentroidYParameter & y
) {
    insert_param(x, params, filter);
    insert_param(y, params, filter);
    return params;     
}

ParamRefs & CentroidParameters::get_parameters(ParamRefs & params, ParamFilter * filter) const {
    return _get_parameters<ParamRefs>(params, filter, this->get_x_param(), this->get_y_param());
}

ParamCRefs & CentroidParameters::get_parameters_const(ParamCRefs & params, ParamFilter * filter) const {
    return _get_parameters<ParamCRefs>(params, filter, this->get_x_param(), this->get_y_param());
}

double CentroidParameters::get_x() const { return _x->get_value(); }
double CentroidParameters::get_y() const { return _y->get_value(); }
std::array<double, 2> CentroidParameters::get_xy() const { return {this->get_x(), this->get_y()}; }

CentroidXParameter & CentroidParameters::get_x_param() const { return *_x; }
CentroidYParameter & CentroidParameters::get_y_param() const { return *_y; }

std::shared_ptr<CentroidXParameter> CentroidParameters::get_x_param_ptr() { return _x; }
std::shared_ptr<CentroidYParameter> CentroidParameters::get_y_param_ptr() { return _y; }

void CentroidParameters::set_x(double x) { _x->set_value(x); }
void CentroidParameters::set_y(double y) { _y->set_value(y); }
void CentroidParameters::set_xy(const std::array<double, 2> & xy) { this->set_x(xy[0]); this->set_y(xy[1]); }

std::string CentroidParameters::str() const {
    return "CentroidParameters(x=" + _x->str() + ", y=" + _y->str() + ")";
}

CentroidParameters::CentroidParameters(
    std::shared_ptr<CentroidXParameter> x,
    std::shared_ptr<CentroidYParameter> y
) :
    CentroidData(),
    _x(x == nullptr ? std::make_shared<CentroidXParameter>() : std::move(x)),
    _y(y == nullptr ? std::make_shared<CentroidYParameter>() : std::move(y))
{}

CentroidParameters::CentroidParameters(double x, double y) : CentroidData(), 
    _x(std::make_shared<CentroidXParameter>(x)),
    _y(std::make_shared<CentroidYParameter>(y))
{}

} // namespace fit
} // namespace gauss2d
