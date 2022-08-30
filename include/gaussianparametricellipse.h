#ifndef GAUSS2D_FIT_GAUSSIANPARAMETRICELLIPSE_H
#define GAUSS2D_FIT_GAUSSIANPARAMETRICELLIPSE_H

#include <memory>

#include "gauss2d/ellipse.h"

#include "parameters.h"
#include "parametricellipse.h"

namespace gauss2d
{
namespace fit
{
class GaussianParametricEllipse : public gauss2d::EllipseData, public ParametricEllipse
{
private:
    std::shared_ptr<SigmaXParameter> _sigma_x;
    std::shared_ptr<SigmaYParameter> _sigma_y;
    std::shared_ptr<RhoParameter> _rho;

public:
    ParamRefs & get_parameters(ParamRefs & params, ParamFilter * filter = nullptr) const override;
    ParamCRefs & get_parameters_const(ParamCRefs & params, ParamFilter * filter = nullptr) const override;

    double get_rho() const override;
    double get_hwhm_x() const override;
    double get_hwhm_y() const override;
    double get_sigma_x() const override;
    double get_sigma_y() const override;
    double get_size_x() const override;
    double get_size_y() const override;

    std::array<double, 3> get_xyr() const override;
    std::array<double, 3> get_hxyr() const override;

    RhoParameter & get_rho_param() const override;
    SigmaXParameter & get_sigma_x_param() const;
    SigmaYParameter & get_sigma_y_param() const;
    SizeXParameter & get_size_x_param() const override;
    SizeYParameter & get_size_y_param() const override;

    std::shared_ptr<RhoParameter> get_rho_param_ptr() override;
    std::shared_ptr<SigmaXParameter> get_sigma_x_param_ptr();
    std::shared_ptr<SigmaYParameter> get_sigma_y_param_ptr();
    std::shared_ptr<SizeXParameter> get_size_x_param_ptr() override;
    std::shared_ptr<SizeYParameter> get_size_y_param_ptr() override;

    void set(double sigma_x, double sigma_y, double rho) override;
    void set_h(double hwhm_x, double hwhm_y, double rho) override;
    void set_rho(double rho) override;
    void set_hwhm_x(double hwhm_x) override;
    void set_hwhm_y(double hwhm_y) override;
    void set_sigma_x(double sigma_x) override;
    void set_sigma_y(double sigma_y) override;
    void set_size_x(double sigma_x) override;
    void set_size_y(double sigma_y) override;
    void set_hxyr(const std::array<double, 3> & hxyr) override;
    void set_xyr(const std::array<double, 3> & xyr) override;

    std::string str() const override;

    GaussianParametricEllipse(
        std::shared_ptr<SigmaXParameter> sigma_x,
        std::shared_ptr<SigmaYParameter> sigma_y,
        std::shared_ptr<RhoParameter> rho=nullptr
    );
    GaussianParametricEllipse(double sigma_x=0, double sigma_y=0, double rho=0);
};
} // namespace fit
} // namespace gauss2d

#endif
