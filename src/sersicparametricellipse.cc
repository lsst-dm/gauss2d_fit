#include "sersicparametricellipse.h"

namespace gauss2d
{
namespace fit
{
template<typename t>
t & _get_parameters(
    t & params, ParamFilter * filter,
    ReffXParameter & x, ReffYParameter & y, RhoParameter & r
) {
    insert_param(x, params, filter);
    insert_param(y, params, filter);
    insert_param(r, params, filter);
    return params;     
}

ParamRefs & SersicParametricEllipse::get_parameters(ParamRefs & params, ParamFilter * filter) const {
    return _get_parameters<ParamRefs>(params, filter, *_size_x, *_size_y, *_rho);
}
ParamCRefs & SersicParametricEllipse::get_parameters_const(ParamCRefs & params, ParamFilter * filter) const {
    return _get_parameters<ParamCRefs>(params, filter, *_size_x, *_size_y, *_rho);
}

double SersicParametricEllipse::get_rho() const { return _rho->get_value(); }
double SersicParametricEllipse::get_size_x() const { return _size_x->get_value(); }
double SersicParametricEllipse::get_size_y() const { return _size_y->get_value(); }
std::array<double, 3> SersicParametricEllipse::get_xyr() const {
    return {this->get_size_x(), this->get_size_y(), this->get_rho()};
}

RhoParameter & SersicParametricEllipse::get_rho_param() const { return *_rho;}
SizeXParameter & SersicParametricEllipse::get_size_x_param() const { return *_size_x; };
SizeYParameter & SersicParametricEllipse::get_size_y_param() const { return *_size_y; };

std::shared_ptr<ReffXParameter> SersicParametricEllipse::get_reff_x_param_ptr() { return _size_x; };
std::shared_ptr<ReffYParameter> SersicParametricEllipse::get_reff_y_param_ptr() { return _size_y; };
std::shared_ptr<RhoParameter> SersicParametricEllipse::get_rho_param_ptr() { return _rho;}
std::shared_ptr<SizeXParameter> SersicParametricEllipse::get_size_x_param_ptr() { return _size_x; };
std::shared_ptr<SizeYParameter> SersicParametricEllipse::get_size_y_param_ptr() { return _size_y; };

inline void SersicParametricEllipse::set_rho(double rho) { _rho->set_value(rho);}
inline void SersicParametricEllipse::set_size_x(double size_x) { _size_x->set_value(size_x); }
inline void SersicParametricEllipse::set_size_y(double size_y) { _size_y->set_value(size_y); }

std::string SersicParametricEllipse::str() const {
    return "SersicParametricEllipse(size_x=" + _size_x->str() + ", size_y=" + _size_y->str()
        + ", rho=" + _rho->str() + ")";
}

SersicParametricEllipse::SersicParametricEllipse(
    std::shared_ptr<ReffXParameter> size_x,
    std::shared_ptr<ReffYParameter> size_y,
    std::shared_ptr<RhoParameter> rho
) :
    _size_x(size_x == nullptr ? std::make_shared<ReffXParameter>(0) : std::move(size_x)),
    _size_y(size_y == nullptr ? std::make_shared<ReffYParameter>(0) : std::move(size_y)),
    _rho(rho == nullptr ? std::make_shared<RhoParameter>(0) : std::move(rho))
{};

SersicParametricEllipse::SersicParametricEllipse(double size_x, double size_y, double rho) :
    _size_x(std::make_shared<ReffXParameter>(size_x)),
    _size_y(std::make_shared<ReffYParameter>(size_y)),
    _rho(std::make_shared<RhoParameter>(rho))
{};

} // namespace fit
} // namespace gauss2d
