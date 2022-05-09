#include "ellipseparameters.h"

namespace gauss2d
{
namespace fit
{
ParamRefs & EllipseParameters::get_parameters(ParamRefs & params, ParamFilter * filter) const {
    insert_param(*_sigma_x, params, filter);
    insert_param(*_sigma_y, params, filter);
    insert_param(*_rho, params, filter);
    return params;
}
ParamCRefs & EllipseParameters::get_parameters_const(ParamCRefs & params, ParamFilter * filter) const {
    insert_param_const(*_sigma_x, params, filter);
    insert_param_const(*_sigma_y, params, filter);
    insert_param_const(*_rho, params, filter);
    return params;
}

double EllipseParameters::get_rho() const { return _rho->get_value(); }
double EllipseParameters::get_sigma_x() const { return _sigma_x->get_value(); }
double EllipseParameters::get_sigma_y() const { return _sigma_y->get_value(); }
std::array<double, 3> EllipseParameters::get_xyr() const {
    return {this->get_sigma_x(), this->get_sigma_y(), this->get_rho()};
}

RhoParameter & EllipseParameters::get_rho_param() const { return *_rho;}
SigmaXParameter & EllipseParameters::get_sigma_x_param() const { return *_sigma_x; }
SigmaYParameter & EllipseParameters::get_sigma_y_param() const { return *_sigma_y; }

std::shared_ptr<RhoParameter> EllipseParameters::get_rho_param_ptr() { return _rho;}
std::shared_ptr<SigmaXParameter> EllipseParameters::get_sigma_x_param_ptr() { return _sigma_x; }
std::shared_ptr<SigmaYParameter> EllipseParameters::get_sigma_y_param_ptr() { return _sigma_y; }

void EllipseParameters::set(double sigma_x, double sigma_y, double rho) {
    set_sigma_x(sigma_x);
    set_sigma_y(sigma_y);
    set_rho(rho);
}
inline void EllipseParameters::set_rho(double rho) { _rho->set_value(rho);}
inline void EllipseParameters::set_sigma_x(double sigma_x) { _sigma_x->set_value(sigma_x); }
inline void EllipseParameters::set_sigma_y(double sigma_y) { _sigma_y->set_value(sigma_y); }
void EllipseParameters::set_xyr(const std::array<double, 3> & xyr) { this->set(xyr[0], xyr[1], xyr[2]); }

std::string EllipseParameters::str() const {
    return "EllipseParameters(sigma_x=" + _sigma_x->str() + ", sigma_y=" + _sigma_y->str()
        + ", rho=" + _rho->str() + ")";
}

EllipseParameters::EllipseParameters(
    std::shared_ptr<SigmaXParameter> sigma_x,
    std::shared_ptr<SigmaYParameter> sigma_y,
    std::shared_ptr<RhoParameter> rho
) :
    _sigma_x(sigma_x == nullptr ? std::make_shared<SigmaXParameter>(0) : std::move(sigma_x)),
    _sigma_y(sigma_y == nullptr ? std::make_shared<SigmaYParameter>(0) : std::move(sigma_y)),
    _rho(rho == nullptr ? std::make_shared<RhoParameter>(0) : std::move(rho))
{};
EllipseParameters::EllipseParameters(double sigma_x, double sigma_y, double rho) :
    _sigma_x(std::make_shared<SigmaXParameter>(sigma_x)),
    _sigma_y(std::make_shared<SigmaYParameter>(sigma_y)),
    _rho(std::make_shared<RhoParameter>(rho))
{};

} // namespace fit
} // namespace gauss2d
