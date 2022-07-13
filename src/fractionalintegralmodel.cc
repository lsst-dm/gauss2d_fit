#include "fractionalintegralmodel.h"

#include <functional>
#include <memory>
#include <stdexcept>
#include <vector>

#include "channel.h"
#include "component.h"
#include "integralmodel.h"
#include "param_defs.h"
#include "param_filter.h"
#include "parameters.h"
#include "util.h"

namespace gauss2d
{
namespace fit
{

ParamRefs & FractionalGaussianIntegral::get_parameters(ParamRefs & params, ParamFilter * filter) const
{
    return _model.get()->get_parameters(params, filter);
}

ParamCRefs & FractionalGaussianIntegral::get_parameters_const(ParamCRefs & params, ParamFilter * filter) const
{
    return _model.get()->get_parameters_const(params, filter);
}

double FractionalGaussianIntegral::get_value() const {
    return _model->get_integral(*_channel);
}

void FractionalGaussianIntegral::set_value(double value) {
    throw std::runtime_error("FractionalGaussianIntegral can't set_value");
}

std::string FractionalGaussianIntegral::str() const {
    return "FractionalGaussianIntegral(channel=" + _channel->str() + ", model="
        + _model->str() + ")";
}

FractionalGaussianIntegral::FractionalGaussianIntegral(
    const std::shared_ptr<Channel> channel,
    std::shared_ptr<FractionalIntegralModel> model
) : _channel(std::move(channel)), _model(std::move(model)) {
    if(_channel == nullptr) throw std::invalid_argument("FractionalGaussianIntegral channel can't be null");
    if(_model == nullptr) throw std::invalid_argument("FractionalGaussianIntegral model can't be null");
    auto channels = _model.get()->get_channels();
    if(channels.find(*channel) == channels.end()) throw std::invalid_argument(
        "FractionalGaussianIntegral channel=" + channel->str() + " not in model.get_channels()="
        + str_iter_refw(channels)
    );
}

std::shared_ptr<ProperFractionParameter> FractionalIntegralModel::at(const Channel & channel) {
    return _data.at(channel);
}
std::shared_ptr<const ProperFractionParameter> FractionalIntegralModel::at(const Channel & channel) const {
    return _data.at(channel);
}

// https://stackoverflow.com/questions/8147027/
// how-do-i-call-stdmake-shared-on-a-class-with-only-protected-or-private-const/
// 8147213#comment58654091_25069711
struct FractionalIntegralModel::Shared_enabler : public FractionalIntegralModel
{
    template <typename... Args>
    Shared_enabler(Args &&... args)
    : FractionalIntegralModel(std::forward<Args>(args)...) {}
};

typename FractionalIntegralModel::Data::iterator FractionalIntegralModel::begin() noexcept {
    return _data.begin();
}
typename FractionalIntegralModel::Data::const_iterator FractionalIntegralModel::cbegin() const noexcept {
    return _data.begin();
}

typename FractionalIntegralModel::Data::iterator FractionalIntegralModel::end() noexcept {
    return _data.end();
}
typename FractionalIntegralModel::Data::const_iterator FractionalIntegralModel::cend() const noexcept {
    return _data.cend();
}

std::set<std::reference_wrapper<const Channel>> FractionalIntegralModel::get_channels() const {
    return map_keys_ref_const(_data);
}

const IntegralModel & FractionalIntegralModel::get_parent_model() const {
    return *_model;
}

inline double _get_integral_remainder(const FractionalIntegralModel * frac, 
    const IntegralModel & model, const Channel & channel) {
    return (frac == nullptr ? model.get_integral(channel) : frac->get_integral_remainder(channel));
}

double FractionalIntegralModel::get_integral(const Channel & channel) const {
    return _data.at(channel)->get_value()*_get_integral_remainder(_parent.get(), *_model, channel);
}

double FractionalIntegralModel::get_integral_remainder(const Channel & channel) const {
    return (1. - _data.at(channel)->get_value())*_get_integral_remainder(_parent.get(), *_model, channel);
}

ParamRefs & FractionalIntegralModel::get_parameters(ParamRefs & params, ParamFilter * filter) const {
    _model->get_parameters(params, filter);
    for(auto & p: _data) params.push_back(*p.second);
    return params;
}

ParamCRefs & FractionalIntegralModel::get_parameters_const(ParamCRefs & params, ParamFilter * filter) const {
    _model->get_parameters_const(params, filter);
    for(const auto & p: _data) params.push_back(*p.second);
    return params;
}

size_t FractionalIntegralModel::size() const { return _data.size(); }

std::string FractionalIntegralModel::str() const {
    std::string s = "FractionalIntegralModel({";
    for(const auto & [channel, integral] : _data) s += channel.get().str() + ": " + integral->str() + ",";
    return s + "})";
}

// Return a pointer to a registered FractionalIntegralModel, if it is one
// TODO: Verify if this is safe, or if dynamic_cast is better
// TODO: Consider a second constructor (but it would likely complicate pybind11'ing)
std::shared_ptr<const FractionalIntegralModel> FractionalIntegralModel::_find_parent(std::shared_ptr<const IntegralModel> model) {
    if(model == nullptr) return nullptr;
    auto parent = std::dynamic_pointer_cast<const FractionalIntegralModel>(model);
    if(parent == nullptr) return nullptr;
    bool found = FractionalIntegralModel::_registry.find(*parent) != _registry.end();
    return found ? parent : nullptr;
}

std::shared_ptr<FractionalIntegralModel> FractionalIntegralModel::make(
    std::optional<const Data> data,
    std::shared_ptr<const IntegralModel> model
) {
    std::shared_ptr<FractionalIntegralModel> ptr = std::make_shared<Shared_enabler>(data, std::move(model));;
    _registry.insert({*ptr, *model});
    FractionalIntegralModel::_registry_rev.insert({*model, ptr});
    return ptr;
}

const std::shared_ptr<const FractionalIntegralModel> FractionalIntegralModel::make_const(
    std::optional<const Data> data,
    std::shared_ptr<const IntegralModel> model
) {
    return make(data, model);
}

FractionalIntegralModel::FractionalIntegralModel(
    std::optional<const Data> data,
    std::shared_ptr<const IntegralModel> model)
    : _model(std::move(model)), _parent(_find_parent(_model))
{
    if(_model == nullptr) throw std::invalid_argument("FractionalIntegralModel model can't be null");
    const auto found = find_model(*_model);
    if(found != nullptr) {
        throw std::invalid_argument("FractionalIntegralModel model=" + _model->str()
            + " already referenced by " + found->str()
        );
    }
    if(data)
    {
        for(const auto & [channel, datum] : *data)
        {
            if(datum == nullptr) throw std::runtime_error("FractionalIntegralModel data["
                + channel.get().str() + "] can't be null");
            _data[channel] = datum;
        }
    } else {
        _data[Channel::NONE()] = std::make_shared<ProperFractionParameter>();
    }
    auto data_keys = map_keys_ref_const(_data);
    auto model_keys = _model.get()->get_channels();
    if(data_keys != model_keys) {
        throw std::invalid_argument(
            "FractionalIntegralModel data channels=" + str_iter_refw(data_keys)
                + " != model.get_channels()=" + str_iter_refw(model_keys)
        );
    }
}
FractionalIntegralModel::~FractionalIntegralModel() {};

} // namespace fit
} // namespace gauss2d
