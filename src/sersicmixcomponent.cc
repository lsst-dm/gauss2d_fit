#include "sersicmixcomponent.h"

#include <stdexcept>
#include <string>

#include "gauss2d/centroid.h"
#include "gauss2d/ellipse.h"
#include "gauss2d/gaussian.h"

#include "channel.h"
#include "component.h"
#include "gaussianmodelintegral.h"
#include "integralmodel.h"
#include "linearsersicmixinterpolator.h"
#include "param_defs.h"
#include "param_filter.h"
#include "sersicmix.h"

namespace gauss2d
{
namespace fit
{

class SersicEllipseData : public EllipseData
{
private:
    const std::shared_ptr<const ReffXParameter> _size_x;
    const std::shared_ptr<const ReffYParameter> _size_y;
    const std::shared_ptr<const RhoParameter> _rho;
    const std::shared_ptr<const SersicMixComponentIndexParameter> _sersicindex;
    unsigned short _index;

    double _get_sizeratio() const { return _sersicindex->get_sizeratio(_index); }

public:
    double get_sigma_x() const override { return _get_sizeratio()*_size_x->get_value(); }
    double get_sigma_y() const override { return _get_sizeratio()*_size_y->get_value(); }
    double get_rho() const override { return _rho->get_value(); }

    void set(double sigma_x, double sigma_y, double rho) override { throw std::runtime_error("Can't set on SersicEllipseData"); }
    void set_h(double hwhm_x, double hwhm_y, double rho) override { throw std::runtime_error("Can't set on SersicEllipseData"); }
    void set_sigma_x(double sigma_x) override { throw std::runtime_error("Can't set on SersicEllipseData"); }
    void set_sigma_y(double sigma_y) override { throw std::runtime_error("Can't set on SersicEllipseData"); }
    void set_rho(double rho) override { throw std::runtime_error("Can't set on SersicEllipseData"); }
    void set_hxyr(const std::array<double, 3> & hxyr) override { throw std::runtime_error("Can't set on SersicEllipseData"); }
    void set_xyr(const std::array<double, 3> & xyr) override { throw std::runtime_error("Can't set on SersicEllipseData"); }

    std::string str() const override {
        return "SersicEllipseData(size_x=" + _size_x->str() + ", size_y=" + _size_y->str() + ", rho=" + _rho->str()
            +  + ")";
    }

    SersicEllipseData(
        const std::shared_ptr<const ReffXParameter> size_x,
        const std::shared_ptr<const ReffYParameter> size_y,
        const std::shared_ptr<const RhoParameter> rho,
        const std::shared_ptr<const SersicMixComponentIndexParameter> sersicindex,
        unsigned short index
    ) :
        _size_x(std::move(size_x)),
        _size_y(std::move(size_y)),
        _rho(std::move(rho)),
        _sersicindex(std::move(sersicindex)),
        _index(index)
    {
        if((_size_x == nullptr) || (_size_y == nullptr) || (_rho == nullptr) || (_sersicindex == nullptr)) {
            throw std::invalid_argument("SersicEllipseData args must not be nullptr");
        }
        if(!(_index < _sersicindex->order)) throw std::invalid_argument("index=" + std::to_string(_index)
            + "!< sersicindex->order=" + std::to_string(_sersicindex->order));
    }
};


class SersicModelIntegral : public GaussianModelIntegral
{
private:
    const std::shared_ptr<const SersicMixComponentIndexParameter> _sersicindex;
    unsigned short _index;

public:
    double get_value() const override {
        return _sersicindex->get_integralratio(_index)*_integralmodel->get_integral(_channel);
    }
    void set_value(double value) override { throw std::runtime_error("Can't set on SersicModelIntegral"); }

    std::string str() const override {
        return "SersicModelIntegral(channel=" + _channel.str() + ", integralmodel=" + _integralmodel->str()
            + ", sersicindex=" + _sersicindex->str() + ", index=" + std::to_string(_index) + ")";
    }

    SersicModelIntegral(
        const Channel & channel, const std::shared_ptr<const IntegralModel> integralmodel,
        const std::shared_ptr<const SersicMixComponentIndexParameter> sersicindex, unsigned short index
    ) : GaussianModelIntegral(channel, integralmodel), _sersicindex(std::move(sersicindex)), _index(index)
    {
        if(_sersicindex == nullptr) throw std::invalid_argument("sersicindex must not be null");
        if(!(_index < _sersicindex->order)) throw std::invalid_argument("index=" + std::to_string(_index)
            + "!< sersicindex->order=" + std::to_string(_sersicindex->order));
    }
    
    ~SersicModelIntegral() {};
};

static const std::shared_ptr<const SersicMixInterpolator> INTERPOLATOR_DEFAULT = std::make_shared<
    const LinearSersicMixInterpolator>(4);

void SersicMixComponentIndexParameter::_set_ratios(double sersicindex) {
    _integralsizes = _interpolator->get_integralsizes(sersicindex);
}

double SersicMixComponentIndexParameter::get_integralratio(unsigned short index) const {
    if(index >= order) throw std::invalid_argument(this->str() + "get_integralratio(index=" + std::to_string(index)
        + " >= max(order=" + std::to_string(order) + "))");
    return _integralsizes[index].integral;
}

double SersicMixComponentIndexParameter::get_sizeratio(unsigned short index) const {
    if(index >= order) throw std::invalid_argument(this->str() + "get_integralratio(index=" + std::to_string(index)
        + " >= max(order=" + std::to_string(order) + "))");
    return _integralsizes[index].sigma;
}

static const std::string limits_sersic_name = std::string(parameters::type_name<
    SersicMixComponentIndexParameter>()) + ".limits_maximal";

static const auto limits_sersic = std::make_shared<const parameters::Limits<double>>(
    0.5, 8.0, parameters::type_name<SersicMixComponentIndexParameter>(), ".limits_maximal"
);

const parameters::Limits<double> & SersicMixComponentIndexParameter::get_limits_maximal() const {
    return *limits_sersic;
}

SersicMixComponentIndexParameter::SersicMixComponentIndexParameter(
    double value,
    std::shared_ptr<const parameters::Limits<double>> limits,
    const std::shared_ptr<const parameters::Transform<double>> transform,
    std::shared_ptr<const parameters::Unit> unit,
    bool fixed,
    std::string label,
    const SetC & inheritors,
    const SetC & modifiers,
    const std::shared_ptr<const SersicMixInterpolator> interpolator
) : SersicIndexParameter(value, nullptr, transform, unit, fixed, label, inheritors, modifiers),
    _interpolator(std::move(interpolator == nullptr ? INTERPOLATOR_DEFAULT : interpolator)),
    order(_interpolator->get_order())
{
    // TODO: determine if this can be avoided
    set_limits(std::move(limits));
    _set_ratios(value);
}

void SersicMixComponentIndexParameter::set_value(double value) {
    SersicIndexParameter::set_value(value);
    _set_ratios(value);
}

void SersicMixComponentIndexParameter::set_value_transformed(double value) {
    SersicIndexParameter::set_value_transformed(value);
    _set_ratios(get_value());
}

// This is the gauss2d convention; see evaluator.h
static const std::array<size_t, 6> IDX_ORDER = {0, 1, 5, 2, 3, 4};

void SersicMixComponent::add_extra_param_map(const Channel & channel, extra_param_map & map, ParameterMap & offsets
    ) const
{
    throw std::runtime_error("not implemented yet");
    map.push_back({0, 0});
}

void SersicMixComponent::add_extra_param_factors(const Channel & channel, extra_param_factors & factors) const
{
    throw std::runtime_error("not implemented yet");
    factors.push_back({0, 0});
}

void SersicMixComponent::add_grad_param_map(const Channel & channel, grad_param_map & map, ParameterMap & offsets
    ) const
{
    throw std::runtime_error("not implemented yet");
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

void SersicMixComponent::add_grad_param_factors(const Channel & channel, grad_param_factors & factors) const
{
    throw std::runtime_error("not implemented yet");
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

std::unique_ptr<const gauss2d::Gaussians> SersicMixComponent::get_gaussians(const Channel & channel) const {
    return std::make_unique<const gauss2d::Gaussians>(_gaussians.at(channel));
}

ParamRefs & SersicMixComponent::get_parameters(ParamRefs & params, ParamFilter * filter) const {
    EllipticalComponent::get_parameters(params, filter);
    insert_param(*_sersicindex, params, filter);
    return params;
}

ParamCRefs & SersicMixComponent::get_parameters_const(ParamCRefs & params, ParamFilter * filter) const {
    EllipticalComponent::get_parameters_const(params, filter);
    insert_param(*_sersicindex, params, filter);
    return params;
}

std::string SersicMixComponent::str() const {
    return "SersicMixComponent(" + EllipticalComponent::str() + ", photo=" + _integralmodel->str() + ")";
}

SersicMixComponent::SersicMixComponent(
    std::shared_ptr<SersicParametricEllipse> ellipse,
    std::shared_ptr<CentroidParameters> centroid,
    std::shared_ptr<IntegralModel> integralmodel,
    std::shared_ptr<SersicMixComponentIndexParameter> sersicindex
) : SersicParametricEllipseHolder(std::move(ellipse)),
    EllipticalComponent(_ellipsedata, centroid, integralmodel),
    _sersicindex(sersicindex != nullptr ? std::move(sersicindex)
        : std::make_shared<SersicMixComponentIndexParameter>())
{
    for(const Channel & channel : _integralmodel->get_channels()) {
        gauss2d::Gaussians::Data gaussians {};
        for(size_t index = 0; index < _sersicindex->order; ++index) {
            auto cen = std::make_shared<Centroid>(this->_centroid);
            auto ell = std::make_shared<Ellipse>(std::make_shared<SersicEllipseData>(
                _ellipsedata->get_reff_x_param_ptr(),
                _ellipsedata->get_reff_y_param_ptr(),
                _ellipsedata->get_rho_param_ptr(),
                _sersicindex,
                index
            ));
            auto integral = std::make_shared<SersicModelIntegral>(
                channel,
                _integralmodel,
                _sersicindex,
                index
            );
            gaussians.emplace_back(std::make_shared<gauss2d::Gaussian>(cen, ell, integral));
        }
        _gaussians[channel] = gaussians;
    };
}

} // namespace fit
} // namespace gauss2d
