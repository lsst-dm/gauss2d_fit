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
    return _data.at(channel)->get_value()*(_get_integral_remainder(_parent.get(), *_model, channel));
}

std::vector<std::pair<ParamBaseCRef, extra_param_factor_values>> 
FractionalIntegralModel::get_integral_derivative_factors(const Channel & channel) const {
    const auto & frac = *(this->_data.at(channel));
    /*
        For a model with no parent:
        weight_comp = weight_total*frac
        dweight_comp/dweight_total = frac

        gauss2d will compute dmodel/dweight_comp
        Fitters will want dmodel/d_frac

        frac = weight_comp/weight_total
        dfrac/dweight_comp = 1/weight_total

        dmodel/d_frac = weight_total*dmodel/dweight_comp

        With a parent:

        parent: weight_comp_parent = weight_total*frac_parent
        child1: weight_comp_child1 = weight_total*(1 - frac_parent)
        ...
        childN: weight_comp_childN = weight_total*(1 - frac_parent)...*(1-frac_childNminus1)

        dmodel/dfrac_parent = weight_total*(1 - frac_parent)*dmodel/dweight_comp
                            = weight_remainder_parent*dmodel/dweight_comp

        i.e. there are N fraction parameters dependent on dmodel/dweight_comp for children.
        (excluding the fixed fraction for an is_final child)
    */
    if(_parent == nullptr) {
        return {{frac, {_model->get_integral(channel), 0., 0.}}};
    }
    auto factors = _parent->get_integral_derivative_factors(channel);
    const double frac_parent = _parent->get_parameter_frac(channel).get_value();
    for(auto & factor : factors) factor.second[0] *= frac_parent;
    if(!is_final()) factors.push_back({frac, {_parent->get_integral_remainder(channel), 0.,}});
    return factors;
}

double FractionalIntegralModel::get_integral_remainder(const Channel & channel) const {
    return (1. - _data.at(channel)->get_value())*_get_integral_remainder(_parent.get(), *_model, channel);
}

ProperFractionParameter & FractionalIntegralModel::get_parameter_frac(const Channel & channel) const {
    return *(_data.at(channel));
}

ParamRefs & FractionalIntegralModel::get_parameters(ParamRefs & params, ParamFilter * filter) const {
    _model->get_parameters(params, filter);
    for(auto & p: _data) insert_param_channel(p.first, *p.second, params, filter);
    return params;
}

ParamCRefs & FractionalIntegralModel::get_parameters_const(ParamCRefs & params, ParamFilter * filter) const {
    _model->get_parameters_const(params, filter);
    for(const auto & p: _data) insert_param_channel(p.first, *p.second, params, filter);
    return params;
}

bool FractionalIntegralModel::is_final() const { return _is_final; }

size_t FractionalIntegralModel::size() const { return _data.size(); }

std::string FractionalIntegralModel::str() const {
    std::string s = "FractionalIntegralModel({";
    for(const auto & [channel, integral] : _data) s += channel.get().str() + ": " + integral->str() + ",";
    return s + "})";
}

// Return a pointer to a registered FractionalIntegralModel, if it is one
// TODO: Verify if this is safe, or if dynamic_cast is better
// TODO: Consider a second constructor (but it would likely complicate pybind11'ing)
std::shared_ptr<const FractionalIntegralModel> FractionalIntegralModel::_find_parent(
    std::shared_ptr<const IntegralModel> model
) {
    if(model == nullptr) return nullptr;
    auto parent = std::dynamic_pointer_cast<const FractionalIntegralModel>(model);
    if(parent == nullptr) return nullptr;
    bool found = FractionalIntegralModel::_registry.find(*parent) != _registry.end();
    return found ? parent : nullptr;
}

std::shared_ptr<FractionalIntegralModel> FractionalIntegralModel::make(
    std::optional<const Data> data,
    std::shared_ptr<const IntegralModel> model,
    bool is_final
) {
    std::shared_ptr<FractionalIntegralModel> ptr = std::make_shared<Shared_enabler>(
        data, std::move(model), is_final);;
    _registry.insert({*ptr, *model});
    FractionalIntegralModel::_registry_rev.insert({*model, ptr});
    return ptr;
}

const std::shared_ptr<const FractionalIntegralModel> FractionalIntegralModel::make_const(
    std::optional<const Data> data,
    std::shared_ptr<const IntegralModel> model,
    bool is_final
) {
    return make(data, model, is_final);
}

FractionalIntegralModel::FractionalIntegralModel(
    std::optional<const Data> data,
    std::shared_ptr<const IntegralModel> model,
    bool is_final)
    : _model(std::move(model)), _parent(_find_parent(_model)), _is_final(is_final)
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
            if(datum == nullptr) {
                throw std::runtime_error("FractionalIntegralModel data["
                    + channel.get().str() + "] can't be null");
            } else if(_is_final) {
                bool is_fixed = datum->get_fixed();
                bool is_one = datum->get_value() == 1;
                std::string errmsg = "";
                if(!is_fixed) errmsg += " is_fixed != true;";
                if(!is_one) errmsg += " get_value()=" + std::to_string(datum->get_value())
                    + "!=1;";
                if(errmsg.size() > 0) {
                    throw std::runtime_error("FractionalIntegralModel data["
                        + channel.get().str() + "] is_final==true but" + errmsg);
                }
            }
            _data[channel] = datum;
        }
    } else {
        auto param = std::make_shared<ProperFractionParameter>();
        param->set_fixed(_is_final);
        _data[Channel::NONE()] = param;
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
