#ifndef GAUSS2D_FIT_SERSICPARAMETRICELLIPSE_H
#define GAUSS2D_FIT_SERSICPARAMETRICELLIPSE_H

#include <memory>

#include "gauss2d/ellipse.h"

#include "parameters.h"
#include "parametricellipse.h"

namespace gauss2d
{
namespace fit
{
class SersicParametricEllipse : public ParametricEllipse
{
private:
    std::shared_ptr<ReffXParameter> _size_x;
    std::shared_ptr<ReffYParameter> _size_y;
    std::shared_ptr<RhoParameter> _rho;

public:
    ParamRefs & get_parameters(ParamRefs & params, ParamFilter * filter = nullptr) const override;
    ParamCRefs & get_parameters_const(ParamCRefs & params, ParamFilter * filter = nullptr) const override;

    double get_rho() const override;
    double get_size_x() const override;
    double get_size_y() const override;
    std::array<double, 3> get_xyr() const override;

    RhoParameter & get_rho_param() const override;
    SizeXParameter & get_size_x_param() const override;
    SizeYParameter & get_size_y_param() const override;

    std::shared_ptr<ReffXParameter> get_reff_x_param_ptr();
    std::shared_ptr<ReffYParameter> get_reff_y_param_ptr();
    std::shared_ptr<RhoParameter> get_rho_param_ptr() override;
    std::shared_ptr<SizeXParameter> get_size_x_param_ptr() override;
    std::shared_ptr<SizeYParameter> get_size_y_param_ptr() override;

    void set_rho(double rho) override;
    void set_size_x(double size_x) override;
    void set_size_y(double size_y) override;

    std::string repr(bool name_keywords = false) const override;
    std::string str() const override;

    SersicParametricEllipse(
        std::shared_ptr<ReffXParameter> size_x,
        std::shared_ptr<ReffYParameter> size_y,
        std::shared_ptr<RhoParameter> rho=nullptr
    );
    SersicParametricEllipse(double size_x=0, double size_y=0, double rho=0);
};
} // namespace fit
} // namespace gauss2d

#endif
