#include "lsst/gauss2d/type_name.h"

#include "lsst/gauss2d/fit/sersicparametricellipse.h"

namespace lsst::gauss2d::fit {
template <typename t>
t& _get_parameters(t& params, ParamFilter* filter, ReffXParameterD& x, ReffYParameterD& y, RhoParameterD& r) {
    insert_param(x, params, filter);
    insert_param(y, params, filter);
    insert_param(r, params, filter);
    return params;
}

SersicParametricEllipse::SersicParametricEllipse(std::shared_ptr<ReffXParameterD> size_x,
                                                 std::shared_ptr<ReffYParameterD> size_y,
                                                 std::shared_ptr<RhoParameterD> rho)
        : _size_x(size_x == nullptr ? std::make_shared<ReffXParameterD>(0) : std::move(size_x)),
          _size_y(size_y == nullptr ? std::make_shared<ReffYParameterD>(0) : std::move(size_y)),
          _rho(rho == nullptr ? std::make_shared<RhoParameterD>(0) : std::move(rho)){};

SersicParametricEllipse::SersicParametricEllipse(double size_x, double size_y, double rho)
        : _size_x(std::make_shared<ReffXParameterD>(size_x)),
          _size_y(std::make_shared<ReffYParameterD>(size_y)),
          _rho(std::make_shared<RhoParameterD>(rho)){};

ParamRefs& SersicParametricEllipse::get_parameters(ParamRefs& params, ParamFilter* filter) const {
    return _get_parameters<ParamRefs>(params, filter, *_size_x, *_size_y, *_rho);
}
ParamCRefs& SersicParametricEllipse::get_parameters_const(ParamCRefs& params, ParamFilter* filter) const {
    return _get_parameters<ParamCRefs>(params, filter, *_size_x, *_size_y, *_rho);
}

double SersicParametricEllipse::get_rho() const { return _rho->get_value(); }
double SersicParametricEllipse::get_size_x() const { return _size_x->get_value(); }
double SersicParametricEllipse::get_size_y() const { return _size_y->get_value(); }
std::array<double, 3> SersicParametricEllipse::get_xyr() const {
    return {this->get_size_x(), this->get_size_y(), this->get_rho()};
}

RhoParameterD& SersicParametricEllipse::get_rho_param() const { return *_rho; }
SizeXParameterD& SersicParametricEllipse::get_size_x_param() const { return *_size_x; };
SizeYParameterD& SersicParametricEllipse::get_size_y_param() const { return *_size_y; };

std::shared_ptr<ReffXParameterD> SersicParametricEllipse::get_reff_x_param_ptr() { return _size_x; };
std::shared_ptr<ReffYParameterD> SersicParametricEllipse::get_reff_y_param_ptr() { return _size_y; };
std::shared_ptr<RhoParameterD> SersicParametricEllipse::get_rho_param_ptr() { return _rho; }
std::shared_ptr<SizeXParameterD> SersicParametricEllipse::get_size_x_param_ptr() { return _size_x; };
std::shared_ptr<SizeYParameterD> SersicParametricEllipse::get_size_y_param_ptr() { return _size_y; };

inline void SersicParametricEllipse::set_rho(double rho) { _rho->set_value(rho); }
inline void SersicParametricEllipse::set_size_x(double size_x) { _size_x->set_value(size_x); }
inline void SersicParametricEllipse::set_size_y(double size_y) { _size_y->set_value(size_y); }

std::string SersicParametricEllipse::repr(bool name_keywords, std::string_view namespace_separator) const {
    return std::string("SersicParametricEllipse(") + (name_keywords ? "size_x=" : "")
           + _size_x->repr(name_keywords, namespace_separator) + ", " + (name_keywords ? "size_y=" : "")
           + _size_y->repr(name_keywords, namespace_separator) + ", " + (name_keywords ? "rho=" : "")
           + _rho->repr(name_keywords, namespace_separator) + ")";
}

std::string SersicParametricEllipse::str() const {
    return "SersicParametricEllipse(size_x=" + _size_x->str() + ", size_y=" + _size_y->str()
           + ", rho=" + _rho->str() + ")";
}
}  // namespace lsst::gauss2d::fit
