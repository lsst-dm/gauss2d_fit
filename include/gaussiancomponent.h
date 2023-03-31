#ifndef GAUSS2D_FIT_GAUSSIANCOMPONENT_H
#define GAUSS2D_FIT_GAUSSIANCOMPONENT_H

#include "channel.h"
#include "ellipticalcomponent.h"
#include "gaussianparametricellipse.h"
#include "integralmodel.h"
#include "linearintegralmodel.h"
#include "param_defs.h"
#include "param_filter.h"
#include <memory>

namespace gauss2d
{
namespace fit
{
// TODO: Revisit the necessity of this class
// Its purpose is to have the GaussianParametricEllipse stored here and initialized first in
// GaussianComponent's constructor
class GaussianParametricEllipseHolder {
public:
    std::shared_ptr<GaussianParametricEllipse> _ellipsedata;

    GaussianParametricEllipseHolder(std::shared_ptr<GaussianParametricEllipse> ellipse = nullptr)
        : _ellipsedata(std::move(ellipse)) {
        if(_ellipsedata == nullptr) _ellipsedata = std::make_shared<GaussianParametricEllipse>();
    }
};

class GaussianComponent : private GaussianParametricEllipseHolder, public EllipticalComponent {
private:
    ParamCRefs _get_parameters_grad(const Channel & channel) const;

public:
    void add_extra_param_map(const Channel & channel, extra_param_map & map_extra,
        const grad_param_map & map_grad, ParameterMap & offsets
    ) const override;
    void add_extra_param_factors(const Channel & channel, extra_param_factors & factors) const override;
    void add_grad_param_map(const Channel & channel, grad_param_map & map, ParameterMap & offsets
        ) const override;
    void add_grad_param_factors(const Channel & channel, grad_param_factors & factor) const override;
    
    std::unique_ptr<const gauss2d::Gaussians> get_gaussians(const Channel & channel) const override;
    size_t get_n_gaussians(const Channel & channel) const override;

    ParamRefs & get_parameters(ParamRefs & params, ParamFilter * filter = nullptr) const override;
    ParamCRefs & get_parameters_const(ParamCRefs & params, ParamFilter * filter = nullptr) const override;
    
    // Useful if you want to make a single-Gaussian, zero-sized component list for the PsfModel
    // when fitting a pre-computed PSF
    static std::vector<std::shared_ptr<Component>> make_uniq_default_gaussians(
        std::vector<double> sizes={2.},
        bool fixed=true
    ) {
        std::vector<std::shared_ptr<Component>> comps = {};
        for(const double size : sizes)
        {
            LinearIntegralModel::Data data = {{
                Channel::NONE(), std::make_shared<IntegralParameter>(1., nullptr, nullptr, nullptr, fixed)
            }};
            comps.emplace_back(std::make_shared<GaussianComponent>(
                std::make_shared<GaussianParametricEllipse>(
                    std::make_shared<SigmaXParameter>(size, nullptr, nullptr, nullptr, fixed),
                    std::make_shared<SigmaYParameter>(size, nullptr, nullptr, nullptr, fixed),
                    std::make_shared<RhoParameter>(0, nullptr, nullptr, nullptr, fixed)
                ),
                std::make_shared<CentroidParameters>(
                    std::make_shared<CentroidXParameter>(0, nullptr, nullptr, nullptr, fixed),
                    std::make_shared<CentroidYParameter>(0, nullptr, nullptr, nullptr, fixed)
                ),
                std::make_shared<LinearIntegralModel>(&data)
            ));
        }
        return comps;
    }

    void set_extra_param_factors(const Channel & channel, extra_param_factors & factors, size_t index) const override;
    void set_grad_param_factors(const Channel & channel, grad_param_factors & factors, size_t index) const override;
    
    std::string repr(bool name_keywords = false) const override;
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
