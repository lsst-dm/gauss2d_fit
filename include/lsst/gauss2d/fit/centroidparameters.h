#ifndef LSST_GAUSS2D_FIT_CENTROIDPARAMETERS_H
#define LSST_GAUSS2D_FIT_CENTROIDPARAMETERS_H

#include <memory>

#include "lsst/gauss2d/centroid.h"

#include "parameters.h"
#include "parametric.h"

namespace lsst::gauss2d::fit {
/**
 * A Centroid with Parameters for x and y
 */
class CentroidParameters : public lsst::gauss2d::CentroidData, public Parametric {
public:
    ParamRefs& get_parameters(ParamRefs& params, ParamFilter* filter = nullptr) const override;
    ParamCRefs& get_parameters_const(ParamCRefs& params, ParamFilter* filter = nullptr) const override;

    double get_x() const override;
    double get_y() const override;
    std::array<double, 2> get_xy() const override;

    /// Get a ref to the x param
    CentroidXParameterD& get_x_param() const;
    /// Get a ref to the y param
    CentroidYParameterD& get_y_param() const;

    std::shared_ptr<CentroidXParameterD> get_x_param_ptr();
    std::shared_ptr<CentroidYParameterD> get_y_param_ptr();

    void set_x(double x) override;
    void set_y(double y) override;
    void set_xy(const std::array<double, 2>& xy) override;

    std::string repr(bool name_keywords = false,
                     std::string_view namespace_separator = Object::CC_NAMESPACE_SEPARATOR) const override;
    std::string str() const override;

    /**
     * Construct a CentroidParameters
     *
     * @param x The x parameter. Default-constructed if null
     * @param y The y parameter. Default-constructed if null
     */
    explicit CentroidParameters(std::shared_ptr<CentroidXParameterD> x = nullptr,
                                std::shared_ptr<CentroidYParameterD> y = nullptr);
    /// Construct a CentroidParameters with default-constructed Parameters set to x/y values
    CentroidParameters(double x, double y);

private:
    std::shared_ptr<CentroidXParameterD> _x;
    std::shared_ptr<CentroidYParameterD> _y;
};
}  // namespace lsst::gauss2d::fit

#endif
