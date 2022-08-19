#ifndef GAUSS2D_FIT_GAUSSIANCOMPONENT_H
#define GAUSS2D_FIT_GAUSSIANCOMPONENT_H

#include "channel.h"
#include "ellipticalcomponent.h"
#include "gaussianparametricellipse.h"
#include "integralmodel.h"
#include "param_defs.h"
#include "param_filter.h"
#include <memory>

namespace gauss2d
{
namespace fit
{
// TODO: Revisit the necessity of this class
class GaussianParametricEllipseHolder {
public:
    std::shared_ptr<GaussianParametricEllipse> _ellipsedata;

    GaussianParametricEllipseHolder(std::shared_ptr<GaussianParametricEllipse> ellipse = nullptr)
        : _ellipsedata(std::move(ellipse)) {
        if(_ellipsedata == nullptr) _ellipsedata = std::make_shared<GaussianParametricEllipse>();
    }
};

class GaussianComponent : private GaussianParametricEllipseHolder, public EllipticalComponent {
public:
    void add_extra_param_map(const Channel & channel, extra_param_map & map, ParameterMap & offsets
        ) const override;
    void add_extra_param_factors(const Channel & channel, extra_param_factors & factors) const override;
    void add_grad_param_map(const Channel & channel, grad_param_map & map, ParameterMap & offsets
        ) const override;
    void add_grad_param_factors(const Channel & channel, grad_param_factors & factor) const override;
    
    std::unique_ptr<const gauss2d::Gaussians> get_gaussians(const Channel & channel) const override;

    ParamRefs & get_parameters(ParamRefs & params, ParamFilter * filter = nullptr) const override;
    ParamCRefs & get_parameters_const(ParamCRefs & params, ParamFilter * filter = nullptr) const override;
    
    std::string str() const override;
    
    GaussianComponent(
        std::shared_ptr<GaussianParametricEllipse> ellipse = nullptr,
        std::shared_ptr<CentroidParameters> centroid = nullptr,
        std::shared_ptr<IntegralModel> integralmodel = nullptr
    );
};
} // namespace fit
} // namespace gauss2d

#endif