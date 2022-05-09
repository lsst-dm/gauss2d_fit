#ifndef GAUSS2D_FIT_CENTROIDPARAMETERS_H
#define GAUSS2D_FIT_CENTROIDPARAMETERS_H

#include <memory>

#include "gauss2d/centroid.h"

#include "parameters.h"
#include "parametric.h"

namespace gauss2d
{
namespace fit
{

class CentroidParameters : public gauss2d::CentroidData, public Parametric
{
private:
    std::shared_ptr<CentroidXParameter> _x;
    std::shared_ptr<CentroidYParameter> _y;

public:
    ParamRefs & get_parameters(ParamRefs & params, ParamFilter * filter = nullptr) const override;
    ParamCRefs & get_parameters_const(ParamCRefs & params, ParamFilter * filter = nullptr) const override;
    
    double get_x() const override;
    double get_y() const override;
    std::array<double, 2> get_xy() const override;

    CentroidXParameter & get_x_param() const;
    CentroidYParameter & get_y_param() const;

    std::shared_ptr<CentroidXParameter> get_x_param_ptr();
    std::shared_ptr<CentroidYParameter> get_y_param_ptr();

    void set_x(double x) override;
    void set_y(double y) override;
    void set_xy(const std::array<double, 2> & xy) override;

    std::string str() const override;

    CentroidParameters(
        std::shared_ptr<CentroidXParameter> x = nullptr,
        std::shared_ptr<CentroidYParameter> y = nullptr
    );

    CentroidParameters(double x, double y);
};
} // namespace fit
} // namespace gauss2d

#endif