#include "gaussiancomponent.h"

#include "gauss2d/gaussian.h"
#include "integralmodel.h"
#include "gaussianparametricellipse.h"
#include "gaussianmodelintegral.h"
#include "param_defs.h"
#include "param_filter.h"
#include <stdexcept>

namespace gauss2d
{
namespace fit
{
// This is the gauss2d convention; see evaluator.h
static const std::array<size_t, 6> IDX_ORDER = {0, 1, 5, 2, 3, 4};

void GaussianComponent::add_extra_param_map(const Channel & channel, extra_param_map & map, ParameterMap & offsets
    ) const
{
    map.push_back({0, 0});
}

void GaussianComponent::add_extra_param_factors(const Channel & channel, extra_param_factors & factors) const
{
    factors.push_back({0, 0, 0});
}

void GaussianComponent::add_grad_param_map(const Channel & channel, grad_param_map & map, ParameterMap & offsets
    ) const
{
    map.push_back({0, 0, 0, 0, 0, 0});
    ParamCRefs params;
    this->get_parameters_const(params);

    auto & values = map.back();

    size_t size_map = offsets.size();
    size_t index_param;
    for(const size_t & order_param : IDX_ORDER) {
        const auto & param = params.at(order_param).get();
        if(!param.get_fixed()) {
            if(offsets.find(param) == offsets.end()) {
                index_param = ++size_map;
                offsets[param] = index_param;
            } else {
                index_param = offsets[param];
            }
            values[order_param] = index_param;
        }
    }
}

void GaussianComponent::add_grad_param_factors(const Channel & channel, grad_param_factors & factors) const
{
    factors.push_back({1, 1, 1, 1, 1, 1});
    ParamCRefs params;
    this->get_parameters_const(params);

    auto & values = factors.back();

    for(const size_t & idx : IDX_ORDER) {
        const auto & param = params.at(idx).get();
        if(param.get_fixed()) {
            values[idx] = 0;
        } else {
            values[idx] *= param.get_transform_derivative();
        }
    }
}

std::unique_ptr<const gauss2d::Gaussians> GaussianComponent::get_gaussians(const Channel & channel) const {
    gauss2d::Gaussians::Data gaussians = {
        std::make_shared<Gaussian>(
            std::make_shared<Centroid>(this->_centroid),
            std::make_shared<Ellipse>(this->_ellipsedata),
            std::make_shared<GaussianModelIntegral>(channel, this->_integralmodel)
        )
    };
    return std::make_unique<const gauss2d::Gaussians>(gaussians);
}

ParamRefs & GaussianComponent::get_parameters(ParamRefs & params, ParamFilter * filter) const {
    EllipticalComponent::get_parameters(params, filter);
    return params;
}

ParamCRefs & GaussianComponent::get_parameters_const(ParamCRefs & params, ParamFilter * filter) const {
    EllipticalComponent::get_parameters_const(params, filter);
    return params;
}

std::string GaussianComponent::str() const {
    return "GaussianComponent(" + EllipticalComponent::str() + ", photo=" + _integralmodel->str() + ")";
}

GaussianComponent::GaussianComponent(
    std::shared_ptr<GaussianParametricEllipse> ellipse,
    std::shared_ptr<CentroidParameters> centroid,
    std::shared_ptr<IntegralModel> integralmodel
) : GaussianParametricEllipseHolder(std::move(ellipse)), EllipticalComponent(_ellipsedata, centroid, integralmodel)
{
}

} // namespace fit
} // namespace gauss2d
