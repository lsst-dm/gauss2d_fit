#ifndef GAUSS2D_FIT_ELLIPSEPARAMETERS_H
#define GAUSS2D_FIT_ELLIPSEPARAMETERS_H

#include <memory>

#include "gauss2d/ellipse.h"
#include "parameters.h"

namespace gauss2d
{
namespace fit
{
class EllipseParameters : public gauss2d::EllipseData
{
private:
    std::shared_ptr<SigmaXParameter> _sigma_x;
    std::shared_ptr<SigmaYParameter> _sigma_y;
    std::shared_ptr<RhoParameter> _rho;

public:
    double get_rho() const { return _rho->get_value(); }
    double get_sigma_x() const { return _sigma_x->get_value(); }
    double get_sigma_y() const { return _sigma_y->get_value(); }
    std::array<double, 3> get_xyr() const;

    RhoParameter & get_rho_param() const { return *_rho;}
    SigmaXParameter & get_sigma_x_param() const { return *_sigma_x; }
    SigmaYParameter & get_sigma_y_param() const { return *_sigma_y; }

    std::shared_ptr<RhoParameter> get_rho_param_ptr() { return _rho;}
    std::shared_ptr<SigmaXParameter> get_sigma_x_param_ptr() { return _sigma_x; }
    std::shared_ptr<SigmaYParameter> get_sigma_y_param_ptr() { return _sigma_y; }

    void set(double sigma_x, double sigma_y, double rho);
    inline void set_rho(double rho) { _rho->set_value(rho);}
    inline void set_sigma_x(double sigma_x) { _sigma_x->set_value(sigma_x); }
    inline void set_sigma_y(double sigma_y) { _sigma_y->set_value(sigma_y); }
    void set_xyr(const std::array<double, 3> & xyr) { this->set(xyr[0], xyr[1], xyr[2]); }

    std::string str() const {
        return "EllipseParameters(sigma_x=" + _rho->str() + ", sigma_y=" + _sigma_x->str() + ", rho=" + _rho->str() + ")";
    };

    EllipseParameters(
        std::shared_ptr<SigmaXParameter> sigma_x,
        std::shared_ptr<SigmaYParameter> sigma_y,
        std::shared_ptr<RhoParameter> rho=nullptr
    ) :
        _sigma_x(sigma_x == nullptr ? std::make_shared<SigmaXParameter>(0) : std::move(sigma_x)),
        _sigma_y(sigma_y == nullptr ? std::make_shared<SigmaYParameter>(0) : std::move(sigma_y)),
        _rho(rho == nullptr ? std::make_shared<RhoParameter>(0) : std::move(rho))
    {};
    EllipseParameters(double sigma_x=0, double sigma_y=0, double rho=0) :
        _sigma_x(std::make_shared<SigmaXParameter>(sigma_x)),
        _sigma_y(std::make_shared<SigmaYParameter>(sigma_y)),
        _rho(std::make_shared<RhoParameter>(rho))
    {};
};
} // namespace fit
} // namespace gauss2d

#endif