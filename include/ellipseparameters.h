#ifndef GAUSS2D_FIT_ELLIPSEPARAMETERS_H
#define GAUSS2D_FIT_ELLIPSEPARAMETERS_H

#include <memory>

#include "gauss2d/ellipse.h"

#include "parameters.h"
#include "parametric.h"

namespace gauss2d
{
namespace fit
{
class EllipseParameters : public gauss2d::EllipseData, public Parametric
{
private:
    std::shared_ptr<SigmaXParameter> _sigma_x;
    std::shared_ptr<SigmaYParameter> _sigma_y;
    std::shared_ptr<RhoParameter> _rho;

public:
    ParamRefs & get_parameters(ParamRefs & params, ParamFilter * filter = nullptr) const override;
    ParamCRefs & get_parameters_const(ParamCRefs & params, ParamFilter * filter = nullptr) const override;

    double get_rho() const;
    double get_sigma_x() const;
    double get_sigma_y() const;
    std::array<double, 3> get_xyr() const;

    RhoParameter & get_rho_param() const;
    SigmaXParameter & get_sigma_x_param() const;
    SigmaYParameter & get_sigma_y_param() const;

    std::shared_ptr<RhoParameter> get_rho_param_ptr();
    std::shared_ptr<SigmaXParameter> get_sigma_x_param_ptr();
    std::shared_ptr<SigmaYParameter> get_sigma_y_param_ptr();

    void set(double sigma_x, double sigma_y, double rho);
    void set_rho(double rho);
    void set_sigma_x(double sigma_x);
    void set_sigma_y(double sigma_y);
    void set_xyr(const std::array<double, 3> & xyr);

    std::string str() const override;

    EllipseParameters(
        std::shared_ptr<SigmaXParameter> sigma_x,
        std::shared_ptr<SigmaYParameter> sigma_y,
        std::shared_ptr<RhoParameter> rho=nullptr
    );
    EllipseParameters(double sigma_x=0, double sigma_y=0, double rho=0);
};
} // namespace fit
} // namespace gauss2d

#endif