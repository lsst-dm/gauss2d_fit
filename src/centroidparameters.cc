#include <memory>

#include "lsst/gauss2d/type_name.h"

#include "lsst/gauss2d/fit/centroidparameters.h"

namespace lsst::gauss2d::fit {
template <typename t>
t& _get_parameters(t& params, ParamFilter* filter, CentroidXParameterD& x, CentroidYParameterD& y) {
    insert_param(x, params, filter);
    insert_param(y, params, filter);
    return params;
}

CentroidParameters::CentroidParameters(std::shared_ptr<CentroidXParameterD> x,
                                       std::shared_ptr<CentroidYParameterD> y)
        : CentroidData(),
          _x(x == nullptr ? std::make_shared<CentroidXParameterD>() : std::move(x)),
          _y(y == nullptr ? std::make_shared<CentroidYParameterD>() : std::move(y)) {}

CentroidParameters::CentroidParameters(double x, double y)
        : CentroidData(),
          _x(std::make_shared<CentroidXParameterD>(x)),
          _y(std::make_shared<CentroidYParameterD>(y)) {}

ParamRefs& CentroidParameters::get_parameters(ParamRefs& params, ParamFilter* filter) const {
    return _get_parameters<ParamRefs>(params, filter, this->get_x_param(), this->get_y_param());
}

ParamCRefs& CentroidParameters::get_parameters_const(ParamCRefs& params, ParamFilter* filter) const {
    return _get_parameters<ParamCRefs>(params, filter, this->get_x_param(), this->get_y_param());
}

double CentroidParameters::get_x() const { return _x->get_value(); }
double CentroidParameters::get_y() const { return _y->get_value(); }
std::array<double, 2> CentroidParameters::get_xy() const { return {this->get_x(), this->get_y()}; }

CentroidXParameterD& CentroidParameters::get_x_param() const { return *_x; }
CentroidYParameterD& CentroidParameters::get_y_param() const { return *_y; }

std::shared_ptr<CentroidXParameterD> CentroidParameters::get_x_param_ptr() { return _x; }
std::shared_ptr<CentroidYParameterD> CentroidParameters::get_y_param_ptr() { return _y; }

void CentroidParameters::set_x(double x) { _x->set_value(x); }
void CentroidParameters::set_y(double y) { _y->set_value(y); }
void CentroidParameters::set_xy(const std::array<double, 2>& xy) {
    this->set_x(xy[0]);
    this->set_y(xy[1]);
}

std::string CentroidParameters::repr(bool name_keywords, std::string_view namespace_separator) const {
    return type_name_str<CentroidParameters>(false, namespace_separator) + "(" + (name_keywords ? "x=" : "")
           + _x->repr(name_keywords, namespace_separator) + ", " + (name_keywords ? "y=" : "")
           + _y->repr(name_keywords, namespace_separator) + ")";
}

std::string CentroidParameters::str() const {
    return type_name_str<CentroidParameters>(true) + "(x=" + _x->str() + ", y=" + _y->str() + ")";
}

}  // namespace lsst::gauss2d::fit
