#ifndef LSST_GAUSS2D_FIT_GAUSSIANPARAMETRICELLIPSE_H
#define LSST_GAUSS2D_FIT_GAUSSIANPARAMETRICELLIPSE_H

#include <memory>

#include "lsst/gauss2d/ellipse.h"

#include "parameters.h"
#include "parametricellipse.h"

namespace lsst::gauss2d::fit {
/**
 * A Parameter-based implementation of lsst::gauss2d::EllipseData and ParametricEllipse.
 */
class GaussianParametricEllipse : public lsst::gauss2d::EllipseData, public ParametricEllipse {
public:
    /**
     * Construct a GaussianParametricEllipse from Parameter instances.
     *
     * @param sigma_x The SigmaXParameter to reference.
     * @param sigma_y The SigmaYParameter to reference.
     * @param rho The RhoParameter to reference. Default-initialized if null.
     */
    explicit GaussianParametricEllipse(std::shared_ptr<SigmaXParameterD> sigma_x,
                                       std::shared_ptr<SigmaYParameterD> sigma_y,
                                       std::shared_ptr<RhoParameterD> rho = nullptr);
    /**
     * Construct a GaussianParametricEllipse from values.
     *
     * @param sigma_x The value of the otherwise default-initialized SigmaXParameter.
     * @param sigma_y The value of the otherwise default-initialized SigmaYParameter.
     * @param rho The value of the otherwise default-initialized RhoParameter.
     */
    explicit GaussianParametricEllipse(double sigma_x = 0, double sigma_y = 0, double rho = 0);

    ParamRefs& get_parameters(ParamRefs& params, ParamFilter* filter = nullptr) const override;
    ParamCRefs& get_parameters_const(ParamCRefs& params, ParamFilter* filter = nullptr) const override;

    double get_rho() const override;
    double get_hwhm_x() const override;
    double get_hwhm_y() const override;
    double get_sigma_x() const override;
    double get_sigma_y() const override;
    double get_size_x() const override;
    double get_size_y() const override;

    std::array<double, 3> get_xyr() const override;
    std::array<double, 3> get_hxyr() const override;

    RhoParameterD& get_rho_param() const override;
    /// Explicit alias for get_size_x_param
    SigmaXParameterD& get_sigma_x_param() const;
    /// Explicit alias for get_size_y_param
    SigmaYParameterD& get_sigma_y_param() const;
    SizeXParameterD& get_size_x_param() const override;
    SizeYParameterD& get_size_y_param() const override;

    std::shared_ptr<RhoParameterD> get_rho_param_ptr() override;
    /// Explicit alias for get_sigma_x_param_ptr
    std::shared_ptr<SigmaXParameterD> get_sigma_x_param_ptr();
    /// Explicit alias for get_sigma_y_param_ptr
    std::shared_ptr<SigmaYParameterD> get_sigma_y_param_ptr();
    std::shared_ptr<SizeXParameterD> get_size_x_param_ptr() override;
    std::shared_ptr<SizeYParameterD> get_size_y_param_ptr() override;

    void set(double sigma_x, double sigma_y, double rho) override;
    void set_h(double hwhm_x, double hwhm_y, double rho) override;
    void set_rho(double rho) override;
    void set_hwhm_x(double hwhm_x) override;
    void set_hwhm_y(double hwhm_y) override;
    void set_sigma_x(double sigma_x) override;
    void set_sigma_y(double sigma_y) override;
    void set_size_x(double sigma_x) override;
    void set_size_y(double sigma_y) override;
    void set_hxyr(const std::array<double, 3>& hxyr) override;
    void set_xyr(const std::array<double, 3>& xyr) override;

    std::string repr(bool name_keywords = false,
                     std::string_view namespace_separator = Object::CC_NAMESPACE_SEPARATOR) const override;
    std::string str() const override;

private:
    std::shared_ptr<SigmaXParameterD> _sigma_x;
    std::shared_ptr<SigmaYParameterD> _sigma_y;
    std::shared_ptr<RhoParameterD> _rho;
};
}  // namespace lsst::gauss2d::fit

#endif
