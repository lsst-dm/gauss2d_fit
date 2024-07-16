#include <stdexcept>

#include "lsst/gauss2d/evaluate.h"
#include "lsst/gauss2d/gaussian.h"
#include "lsst/gauss2d/type_name.h"

#include "lsst/gauss2d/fit/gaussiancomponent.h"
#include "lsst/gauss2d/fit/integralmodel.h"
#include "lsst/gauss2d/fit/gaussianparametricellipse.h"
#include "lsst/gauss2d/fit/gaussianmodelintegral.h"
#include "lsst/gauss2d/fit/param_defs.h"
#include "lsst/gauss2d/fit/param_filter.h"

namespace lsst::gauss2d::fit {
// This is the gauss2d convention; see evaluator.h
static const std::array<size_t, N_PARAMS_GAUSS2D> IDX_ORDER = {0, 1, 3, 4, 5, 2};
static const size_t N_PARAMS_INTEGRAL_MAX = 2;
static const size_t N_PARAMS_EXTRA_INTEGRAL_MAX = N_PARAMS_INTEGRAL_MAX - 1;

GaussianComponent::GaussianComponent(std::shared_ptr<GaussianParametricEllipse> ellipse,
                                     std::shared_ptr<CentroidParameters> centroid,
                                     std::shared_ptr<IntegralModel> integralmodel)
        : GaussianParametricEllipseHolder(std::move(ellipse)),
          EllipticalComponent(_ellipsedata, centroid, integralmodel) {}

ParamCRefs GaussianComponent::_get_parameters_grad(const Channel& channel) const {
    ParamCRefs params;
    ParamFilter filter{};
    filter.channel = channel;
    _centroid->get_parameters_const(params, &filter);
    _ellipse->get_parameters_const(params, &filter);

    const size_t n_params_nonint = params.size();
    size_t n_params_expect = N_PARAMS_GAUSS2D - 1;

    if (n_params_nonint != n_params_expect) {
        throw std::runtime_error(this->str() + " (centroid + ellipse).get_parameters_const(channel="
                                 + channel.str() + ").size())=" + std::to_string(params.size())
                                 + "!=(N_PARAMS_GAUSS2D - 1)=" + std::to_string(n_params_expect));
    }

    ParamCRefs params_int;
    this->_integralmodel->get_parameters_const(params_int, &filter);

    size_t n_params_int = params_int.size();
    if (!(n_params_int <= N_PARAMS_INTEGRAL_MAX)) {
        throw std::runtime_error(this->str() + "->_integralmodel->get_parameters_const(channel="
                                 + channel.str() + ").size())=" + std::to_string(params_int.size())
                                 + "!<=" + std::to_string(N_PARAMS_INTEGRAL_MAX));
    }
    for (auto param : params_int) params.push_back(param);
    return params;
}

void GaussianComponent::add_extra_param_map(const Channel& channel, ExtraParamMap& map_extra,
                                            const GradParamMap& map_grad, ParameterMap& offsets) const {
    auto factors_int = _integralmodel->get_integral_derivative_factors(channel);
    const size_t n_factors = factors_int.size();
    if (!(n_factors <= N_PARAMS_EXTRA_INTEGRAL_MAX)) {
        throw std::runtime_error(_integralmodel->str() + ".get_integral_derivative_factors(" + channel.str()
                                 + ").size()=" + std::to_string(n_factors)
                                 + "!<=" + std::to_string(N_PARAMS_EXTRA_INTEGRAL_MAX));
    }
    size_t idx = 0;
    size_t offset = 0;

    for (const auto& factor : factors_int) {
        const auto& param = factor.first.get();
        if (!param.get_fixed()) {
            auto found = offsets.find(param);
            if (found == offsets.end()) {
                throw std::runtime_error("param=" + param.str() + " not found in offsets;"
                    " was add_grad_param_map called first?");
            }
            offset = (*found).second;
            idx = map_grad.size() - 1;
        }
    }
    map_extra.push_back({idx, offset});
}

void GaussianComponent::add_extra_param_factors(const Channel& channel, ExtraParamFactors& factors) const {
    factors.push_back({0, 0, 0});
}

void GaussianComponent::add_grad_param_map(const Channel& channel, GradParamMap& map,
                                           ParameterMap& offsets) const {
    auto params = nonconsecutive_unique(this->_get_parameters_grad(channel));
    map.push_back({0, 0, 0, 0, 0, 0});
    auto& values = map.back();

    size_t size_map = offsets.size();
    size_t index_param_map;
    for (size_t idx_param = 0; idx_param < N_PARAMS_GAUSS2D; ++idx_param) {
        const size_t& order_param = IDX_ORDER[idx_param];
        // The parameters must be in the same order as returned by get_parameters(_const)
        const auto& param = params.at(idx_param).get();
        if (!param.get_fixed()) {
            if (offsets.find(param) == offsets.end()) {
                index_param_map = ++size_map;
                offsets[param] = index_param_map;
            } else {
                index_param_map = offsets[param];
            }
            values[order_param] = index_param_map;
        }
    }
    const size_t n_params = params.size();
    // Add offset for any extra integral parameters
    for (size_t idx_param = N_PARAMS_GAUSS2D; idx_param < n_params; ++idx_param) {
        const auto& param = params.at(idx_param).get();
        if (!param.get_fixed() && (offsets.find(param) == offsets.end())) {
            offsets[param] = ++size_map;
        }
    }
}

void GaussianComponent::add_grad_param_factors(const Channel& channel, GradParamFactors& factors) const {
    factors.push_back({1, 1, 1, 1, 1, 1});
}

std::unique_ptr<const lsst::gauss2d::Gaussians> GaussianComponent::get_gaussians(
        const Channel& channel) const {
    lsst::gauss2d::Gaussians::Data gaussians = {std::make_shared<Gaussian>(
            std::make_shared<Centroid>(this->_centroid), std::make_shared<Ellipse>(this->_ellipsedata),
            std::make_shared<GaussianModelIntegral>(channel, this->_integralmodel))};
    return std::make_unique<const lsst::gauss2d::Gaussians>(gaussians);
}

size_t GaussianComponent::get_n_gaussians(const Channel& channel) const { return 1; };

ParamRefs& GaussianComponent::get_parameters(ParamRefs& params, ParamFilter* filter) const {
    EllipticalComponent::get_parameters(params, filter);
    return params;
}

ParamCRefs& GaussianComponent::get_parameters_const(ParamCRefs& params, ParamFilter* filter) const {
    EllipticalComponent::get_parameters_const(params, filter);
    return params;
}

void GaussianComponent::set_extra_param_factors(const Channel& channel, ExtraParamFactors& factors,
                                                size_t index) const {
    auto factors_int = _integralmodel->get_integral_derivative_factors(channel);
    const size_t n_factors = factors_int.size();
    if (!(n_factors <= 1)) {
        throw std::runtime_error(this->str() + ".set_extra_param_factors(" + channel.str()
                                 + ", index=" + std::to_string(index) + " has " + std::to_string(n_factors)
                                 + " parameter factors" + "; only a maximum of one is supported.");
    }
    auto& row = factors[index];

    if (n_factors > 0) {
        row[0] = factors_int[0].second[0] / factors_int[0].first.get().get_transform_derivative();
    } else {
        row[0] = 0;
    }
    row[1] = 0;
    row[2] = 0;
}

void GaussianComponent::set_grad_param_factors(const Channel& channel, GradParamFactors& factors,
                                               size_t index) const {
    auto params = this->_get_parameters_grad(channel);
    auto& values = factors.at(index);

    for (size_t idx_param = 0; idx_param < N_PARAMS_GAUSS2D; ++idx_param) {
        const size_t& order_param = IDX_ORDER[idx_param];
        // The parameters must be in the same order as returned by get_parameters(_const)
        const auto& param = params[idx_param].get();
        if (param.get_fixed()) {
            values[order_param] = 0;
        } else {
            auto deriv = param.get_transform_derivative();
            if (deriv == 0) {
                throw std::runtime_error("Param[idx=" + std::to_string(idx_param) + "]=" + param.str()
                                         + " get_transform_derivative=0 (will result in divide by 0)");
            }
            values[order_param] = 1.0 / deriv;
        }
    }
}

std::string GaussianComponent::repr(bool name_keywords, std::string_view namespace_separator) const {
    return type_name_str<GaussianComponent>(false, namespace_separator) + "("
           + EllipticalComponent::repr(name_keywords, namespace_separator) + ")";
}

std::string GaussianComponent::str() const {
    return type_name_str<GaussianComponent>(true) + "(" + EllipticalComponent::str() + ")";
}

}  // namespace lsst::gauss2d::fit
