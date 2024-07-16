#ifndef LSST_GAUSS2D_FIT_SERSICPARAMETRICELLIPSE_H
#define LSST_GAUSS2D_FIT_SERSICPARAMETRICELLIPSE_H

#include <memory>

#include "lsst/gauss2d/ellipse.h"

#include "parameters.h"
#include "parametricellipse.h"

namespace lsst::gauss2d::fit {
/**
 * A ParametericEllipse with effective radius Parameters.
 */
class SersicParametricEllipse : public ParametricEllipse {
public:
    /**
     * Construct a SersicParametricEllipse with existing Parameter instances.
     *
     * @param size_x The x-axis effective radius parameter.
     * @param size_y The y-axis effective radius parameter.
     * @param rho The correlation (rho) parameter.
     */
    explicit SersicParametricEllipse(std::shared_ptr<ReffXParameterD> size_x,
                                     std::shared_ptr<ReffYParameterD> size_y,
                                     std::shared_ptr<RhoParameterD> rho = nullptr);
    explicit SersicParametricEllipse(double size_x = 0, double size_y = 0, double rho = 0);

    ParamRefs& get_parameters(ParamRefs& params, ParamFilter* filter = nullptr) const override;
    ParamCRefs& get_parameters_const(ParamCRefs& params, ParamFilter* filter = nullptr) const override;

    double get_rho() const override;
    double get_size_x() const override;
    double get_size_y() const override;
    std::array<double, 3> get_xyr() const override;

    RhoParameterD& get_rho_param() const override;
    SizeXParameterD& get_size_x_param() const override;
    SizeYParameterD& get_size_y_param() const override;

    std::shared_ptr<ReffXParameterD> get_reff_x_param_ptr();
    std::shared_ptr<ReffYParameterD> get_reff_y_param_ptr();
    std::shared_ptr<RhoParameterD> get_rho_param_ptr() override;
    std::shared_ptr<SizeXParameterD> get_size_x_param_ptr() override;
    std::shared_ptr<SizeYParameterD> get_size_y_param_ptr() override;

    void set_rho(double rho) override;
    void set_size_x(double size_x) override;
    void set_size_y(double size_y) override;

    std::string repr(bool name_keywords = false,
                     std::string_view namespace_separator = Object::CC_NAMESPACE_SEPARATOR) const override;
    std::string str() const override;

private:
    std::shared_ptr<ReffXParameterD> _size_x;
    std::shared_ptr<ReffYParameterD> _size_y;
    std::shared_ptr<RhoParameterD> _rho;
};
}  // namespace lsst::gauss2d::fit

#endif
