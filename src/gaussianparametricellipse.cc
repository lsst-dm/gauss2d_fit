#include "lsst/gauss2d/type_name.h"

#include "lsst/gauss2d/ellipse.h"
#include "lsst/gauss2d/fit/gaussianparametricellipse.h"
#include "lsst/gauss2d/fit/parameters.h"

namespace lsst::gauss2d::fit {
template <typename t>
t& _get_parameters(t& params, ParamFilter* filter, SigmaXParameterD& x, SigmaYParameterD& y,
                   RhoParameterD& r) {
    insert_param(x, params, filter);
    insert_param(y, params, filter);
    insert_param(r, params, filter);
    return params;
}

GaussianParametricEllipse::GaussianParametricEllipse(std::shared_ptr<SigmaXParameterD> sigma_x,
                                                     std::shared_ptr<SigmaYParameterD> sigma_y,
                                                     std::shared_ptr<RhoParameterD> rho)
        : _sigma_x(sigma_x == nullptr ? std::make_shared<SigmaXParameterD>(0) : std::move(sigma_x)),
          _sigma_y(sigma_y == nullptr ? std::make_shared<SigmaYParameterD>(0) : std::move(sigma_y)),
          _rho(rho == nullptr ? std::make_shared<RhoParameterD>(0) : std::move(rho)){};
GaussianParametricEllipse::GaussianParametricEllipse(double sigma_x, double sigma_y, double rho)
        : _sigma_x(std::make_shared<SigmaXParameterD>(sigma_x)),
          _sigma_y(std::make_shared<SigmaYParameterD>(sigma_y)),
          _rho(std::make_shared<RhoParameterD>(rho)){};

ParamRefs& GaussianParametricEllipse::get_parameters(ParamRefs& params, ParamFilter* filter) const {
    return _get_parameters<ParamRefs>(params, filter, *_sigma_x, *_sigma_y, *_rho);
}

ParamCRefs& GaussianParametricEllipse::get_parameters_const(ParamCRefs& params, ParamFilter* filter) const {
    return _get_parameters<ParamCRefs>(params, filter, *_sigma_x, *_sigma_y, *_rho);
}

double GaussianParametricEllipse::get_hwhm_x() const {
    return lsst::gauss2d::M_SIGMA_HWHM * _sigma_x->get_value();
}
double GaussianParametricEllipse::get_hwhm_y() const {
    return lsst::gauss2d::M_SIGMA_HWHM * _sigma_y->get_value();
}
double GaussianParametricEllipse::get_rho() const { return _rho->get_value(); }
double GaussianParametricEllipse::get_sigma_x() const { return _sigma_x->get_value(); }
double GaussianParametricEllipse::get_sigma_y() const { return _sigma_y->get_value(); }
double GaussianParametricEllipse::get_size_x() const { return _sigma_x->get_value(); }
double GaussianParametricEllipse::get_size_y() const { return _sigma_y->get_value(); }
std::array<double, 3> GaussianParametricEllipse::get_hxyr() const {
    return {this->get_hwhm_x(), this->get_hwhm_y(), this->get_rho()};
}

std::array<double, 3> GaussianParametricEllipse::get_xyr() const {
    return {this->get_sigma_x(), this->get_sigma_y(), this->get_rho()};
}

RhoParameterD& GaussianParametricEllipse::get_rho_param() const { return *_rho; }
SigmaXParameterD& GaussianParametricEllipse::get_sigma_x_param() const { return *_sigma_x; }
SigmaYParameterD& GaussianParametricEllipse::get_sigma_y_param() const { return *_sigma_y; }
SizeXParameterD& GaussianParametricEllipse::get_size_x_param() const { return *_sigma_x; };
SizeYParameterD& GaussianParametricEllipse::get_size_y_param() const { return *_sigma_y; };

std::shared_ptr<RhoParameterD> GaussianParametricEllipse::get_rho_param_ptr() { return _rho; }
std::shared_ptr<SigmaXParameterD> GaussianParametricEllipse::get_sigma_x_param_ptr() { return _sigma_x; }
std::shared_ptr<SigmaYParameterD> GaussianParametricEllipse::get_sigma_y_param_ptr() { return _sigma_y; }
std::shared_ptr<SizeXParameterD> GaussianParametricEllipse::get_size_x_param_ptr() { return _sigma_x; };
std::shared_ptr<SizeYParameterD> GaussianParametricEllipse::get_size_y_param_ptr() { return _sigma_y; };

void GaussianParametricEllipse::set(double sigma_x, double sigma_y, double rho) {
    set_sigma_x(sigma_x);
    set_sigma_y(sigma_y);
    set_rho(rho);
}
void GaussianParametricEllipse::set_h(double hwhm_x, double hwhm_y, double rho) {
    set_hwhm_x(hwhm_x);
    set_hwhm_y(hwhm_y);
    set_rho(rho);
}
inline void GaussianParametricEllipse::set_hwhm_x(double hwhm_x) {
    _sigma_x->set_value(lsst::gauss2d::M_HWHM_SIGMA * hwhm_x);
}
inline void GaussianParametricEllipse::set_hwhm_y(double hwhm_y) {
    _sigma_y->set_value(lsst::gauss2d::M_HWHM_SIGMA * hwhm_y);
}
inline void GaussianParametricEllipse::set_rho(double rho) { _rho->set_value(rho); }
inline void GaussianParametricEllipse::set_sigma_x(double sigma_x) { _sigma_x->set_value(sigma_x); }
inline void GaussianParametricEllipse::set_sigma_y(double sigma_y) { _sigma_y->set_value(sigma_y); }
inline void GaussianParametricEllipse::set_size_x(double size_x) { _sigma_x->set_value(size_x); }
inline void GaussianParametricEllipse::set_size_y(double size_y) { _sigma_y->set_value(size_y); }
void GaussianParametricEllipse::set_hxyr(const std::array<double, 3>& hxyr) {
    this->set_h(hxyr[0], hxyr[1], hxyr[2]);
}
void GaussianParametricEllipse::set_xyr(const std::array<double, 3>& xyr) {
    this->set(xyr[0], xyr[1], xyr[2]);
}

std::string GaussianParametricEllipse::repr(bool name_keywords, std::string_view namespace_separator) const {
    return type_name_str<GaussianParametricEllipse>(false, namespace_separator) + "("
           + (name_keywords ? "sigma_x=" : "") + _sigma_x->repr(name_keywords, namespace_separator) + ", "
           + (name_keywords ? "sigma_y=" : "") + _sigma_y->repr(name_keywords, namespace_separator) + ", "
           + (name_keywords ? "rho=" : "") + _rho->repr(name_keywords, namespace_separator) + ")";
}

std::string GaussianParametricEllipse::str() const {
    return type_name_str<GaussianParametricEllipse>(true) + "(sigma_x=" + _sigma_x->str()
           + ", sigma_y=" + _sigma_y->str() + ", rho=" + _rho->str() + ")";
}

}  // namespace lsst::gauss2d::fit
