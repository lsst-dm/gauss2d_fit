#include "gaussiancomponent.h"

#include "gauss2d/evaluate.h"
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
static const std::array<size_t, N_PARAMS_GAUSS2D> IDX_ORDER = {0, 1, 3, 4, 5, 2};

void GaussianComponent::add_extra_param_map(
    const Channel & channel, extra_param_map & map_extra, const grad_param_map & map_grad, ParameterMap & offsets
) const
{
    map_extra.push_back({0, 0});
}

void GaussianComponent::add_extra_param_factors(const Channel & channel, extra_param_factors & factors) const
{
    factors.push_back({0, 0, 0});
}

void GaussianComponent::add_grad_param_map(
    const Channel & channel, grad_param_map & map, ParameterMap & offsets
) const
{
    ParamCRefs params;
    ParamFilter filter{};
    filter.channel = channel;
    _centroid->get_parameters_const(params, &filter);
    _ellipse->get_parameters_const(params, &filter);

    const size_t n_params_nonint = params.size();
    size_t n_params_expect = N_PARAMS_GAUSS2D - 1;

    if(n_params_nonint != n_params_expect) {
        throw std::runtime_error(
        this->str() + " (centroid + ellipse).get_parameters_const(channel=" + channel.str() + ").size())="
            + std::to_string(params.size()) + "!=(N_PARAMS_GAUSS2D - 1)="
            + std::to_string(n_params_expect)
        );
    }
    
    ParamCRefs params_int;
    filter.nonlinear = false;
    this->_integralmodel->get_parameters_const(params_int, &filter);

    size_t n_params_int = params_int.size();
    n_params_expect = 1;
    if(n_params_int != n_params_expect) {
        throw std::runtime_error(
        this->str() + "->_integralmodel->get_parameters_const(channel=" + channel.str() + ").size())="
            + std::to_string(params_int.size()) + "!=" + std::to_string(n_params_expect)
        );
    }

    params.push_back(params_int[0]);

    map.push_back({0, 0, 0, 0, 0, 0});
    auto & values = map.back();

    size_t size_map = offsets.size();
    size_t index_param_map;
    for(size_t idx_param = 0; idx_param < N_PARAMS_GAUSS2D; ++idx_param) {
        const size_t & order_param = IDX_ORDER[idx_param];
        // The parameters must be in the same order as returned by get_parameters(_const)
        const auto & param = params.at(idx_param).get();
        if(!param.get_fixed()) {
            if(offsets.find(param) == offsets.end()) {
                index_param_map = ++size_map;
                offsets[param] = index_param_map;
            } else {
                index_param_map = offsets[param];
            }
            values[order_param] = index_param_map;
        }
    }
}

void GaussianComponent::add_grad_param_factors(const Channel & channel, grad_param_factors & factors) const
{
    factors.push_back({1, 1, 1, 1, 1, 1});
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

size_t GaussianComponent::get_n_gaussians(const Channel & channel) const { return 1; };

ParamRefs & GaussianComponent::get_parameters(ParamRefs & params, ParamFilter * filter) const {
    EllipticalComponent::get_parameters(params, filter);
    return params;
}

ParamCRefs & GaussianComponent::get_parameters_const(ParamCRefs & params, ParamFilter * filter) const {
    EllipticalComponent::get_parameters_const(params, filter);
    return params;
}

void GaussianComponent::set_extra_param_factors(
    const Channel & channel, extra_param_factors & factors, size_t index
) const {
    auto factors_int = _integralmodel->get_integral_derivative_factors(channel);
    const size_t n_factors = factors_int.size();
    if(n_factors > 1) {
        throw std::runtime_error(this->str() + ".set_extra_param_factors(" + channel.str() + ", index="
            + std::to_string(index) + " has " + std::to_string(n_factors) + " parameter factors"
            + "; only a maximum of one is supported."
        );
    }
    auto & row = factors[index];
    
    if(n_factors == 0) {
        row[0] = 0;
        row[1] = 0;
        row[2] = 0;
    } else {

    }
}

void GaussianComponent::set_grad_param_factors(
    const Channel & channel, grad_param_factors & factors, size_t index
) const {
    ParamCRefs params;
    ParamFilter filter{};
    filter.channel = channel;
    _centroid->get_parameters_const(params, &filter);
    _ellipse->get_parameters_const(params, &filter);
    filter.nonlinear = false;
    this->_integralmodel->get_parameters_const(params, &filter);

    if(params.size() != N_PARAMS_GAUSS2D) {
        throw std::runtime_error(
            this->str() + ".get_parameters_const(" + channel.str() + ").size()="
                + std::to_string(params.size()) + "!=N_PARAMS_GAUSS2D=" + std::to_string(N_PARAMS_GAUSS2D)
                + " (excluding nonlinear integralmodel params)"
        );
    }

    auto & values = factors.at(index);

    for(size_t idx_param = 0; idx_param < N_PARAMS_GAUSS2D; ++idx_param) {
        const size_t & order_param = IDX_ORDER[idx_param];
        // The parameters must be in the same order as returned by get_parameters(_const)
        const auto & param = params.at(idx_param).get();
        if(param.get_fixed()) {
            values[order_param] = 0;
        } else {
            const auto deriv = param.get_transform_derivative();
            if(deriv == 0) throw std::runtime_error("Param[idx=" + std::to_string(idx_param) + "="
                + param.str() + " get_transform_derivative=0 (will result in divide by 0)");
            values[order_param] /= deriv;
        }
    }
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
