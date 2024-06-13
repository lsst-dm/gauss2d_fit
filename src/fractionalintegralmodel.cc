#include <functional>
#include <memory>
#include <stdexcept>
#include <vector>

#include "lsst/gauss2d/to_string.h"
#include "lsst/gauss2d/type_name.h"

#include "lsst/gauss2d/fit/channel.h"
#include "lsst/gauss2d/fit/component.h"
#include "lsst/gauss2d/fit/fractionalintegralmodel.h"
#include "lsst/gauss2d/fit/integralmodel.h"
#include "lsst/gauss2d/fit/param_defs.h"
#include "lsst/gauss2d/fit/param_filter.h"
#include "lsst/gauss2d/fit/parameters.h"
#include "lsst/gauss2d/fit/util.h"

namespace lsst::gauss2d::fit {

FractionalIntegralModel::FractionalIntegralModel(std::optional<const Data> data,
                                                 std::shared_ptr<const IntegralModel> model, bool is_final)
        : _model(std::move(model)), _parent(_find_parent(_model)), _is_final(is_final) {
    if (_model == nullptr) throw std::invalid_argument("FractionalIntegralModel model can't be null");

    const auto found = find_model(*_model);
    if (found != nullptr) {
        throw std::invalid_argument("FractionalIntegralModel model=" + _model->str()
                                    + " already referenced by " + found->str());
    }
    if (data) {
        size_t idx = 0;
        for (const auto &datum : *data) {
            const auto &channel = datum.first;
            const auto &param = datum.second;
            if (param == nullptr) {
                throw std::runtime_error("FractionalIntegralModel data[" + channel.get().str()
                                         + "] can't be null");
            } else if (_is_final) {
                bool is_fixed = param->get_fixed();
                bool is_one = param->get_value() == 1;
                std::string errmsg = "";
                if (!is_fixed) errmsg += " is_fixed != true;";
                if (!is_one) errmsg += " get_value()=" + std::to_string(param->get_value()) + "!=1;";
                if (errmsg.size() > 0) {
                    throw std::invalid_argument("FractionalIntegralModel data[" + std::to_string(idx)
                                                + "] is_final==true but param for " + channel.get().str()
                                                + errmsg);
                }
            }
            if (_map.find(channel) != _map.end()) {
                throw std::runtime_error("FractionalIntegralModel data[" + std::to_string(idx)
                                         + "] channel=" + channel.get().str() + " duplicated");
            }
            _data.emplace_back(datum.first, datum.second);
            _map.insert(_data.back());
            idx++;
        }
    } else {
        auto param = std::make_shared<ProperFractionParameterD>();
        param->set_fixed(_is_final);
        _data.emplace_back(Channel::NONE(), param);
        _map.insert(_data.back());
    }
    auto data_keys = this->get_channels();
    auto model_keys = _model.get()->get_channels();
    if (data_keys != model_keys) {
        throw std::invalid_argument("FractionalIntegralModel data channels=" + str_iter_ref<true>(data_keys)
                                    + " != model.get_channels()=" + str_iter_ref<true>(model_keys));
    }
}
FractionalIntegralModel::~FractionalIntegralModel(){};

std::shared_ptr<ProperFractionParameterD> FractionalIntegralModel::at(const Channel &channel) {
    return _map.at(channel);
}
std::shared_ptr<const ProperFractionParameterD> FractionalIntegralModel::at(const Channel &channel) const {
    return _map.at(channel);
}

// https://stackoverflow.com/questions/8147027/
// how-do-i-call-stdmake-shared-on-a-class-with-only-protected-or-private-const/
// 8147213#comment58654091_25069711
struct FractionalIntegralModel::Shared_enabler : public FractionalIntegralModel {
    template <typename... Args>
    Shared_enabler(Args &&...args) : FractionalIntegralModel(std::forward<Args>(args)...) {}
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

std::vector<std::reference_wrapper<const Channel>> FractionalIntegralModel::get_channels() const {
    std::vector<std::reference_wrapper<const Channel>> rval = {};
    for (auto &datum : _data) rval.emplace_back(datum.first);
    return rval;
}

const IntegralModel &FractionalIntegralModel::get_parent_model() const { return *_model; }

inline double _get_integral_remainder(const FractionalIntegralModel *frac, const IntegralModel &model,
                                      const Channel &channel) {
    return (frac == nullptr ? model.get_integral(channel) : frac->get_integral_remainder(channel));
}

double FractionalIntegralModel::get_integral(const Channel &channel) const {
    return _map.at(channel)->get_value() * (_get_integral_remainder(_parent.get(), *_model, channel));
}

std::vector<std::pair<ParamBaseCRef, ExtraParamFactorValues>>
FractionalIntegralModel::get_integral_derivative_factors(const Channel &channel) const {
    const auto &frac = *(this->_map.at(channel));
    /*
        For a model with no parent:
        gauss2d will evaluate dmodel/dweight_comp
        Fitters will want dmodel/dfrac

        frac = weight_comp/weight_total
        dfrac/dweight_comp = 1/weight_total
        dweight_comp = dfrac*weight_total

        dmodel/d_frac = weight_total*dmodel/dweight_comp

        With a parent:

        parent: weight_comp_parent = weight_total*frac_parent
        child1: weight_comp_child1 = weight_total*(1 - frac_parent)
        ...
        childN: weight_comp_childN = weight_total*(1 - frac_parent)...*(1-frac_childNminus1)

        dmodel/dfrac_parent = weight_total*(1 - frac_parent)*dmodel/dweight_comp
                            = -weight_total*dmodel/dweight_comp

        i.e. there are N fraction parameters dependent on dmodel/dweight_comp for children.
        (excluding the fixed fraction for an is_final child)
    */
    if (_parent == nullptr) {
        return {{frac, {_model->get_integral(channel), 0., 0.}}};
    }
    auto factors = _parent->get_integral_derivative_factors(channel);
    if (is_final()) {
        /*
            The last component's fraction is fixed, but it has an integral of
            (1 - frac_previous)*integral_remaining, the derivative of which
            w.r.t. frac_revious is -1*integral_remaining.
         */
        factors.back().second[0] *= -1.;
    } else {
        // TODO: Check this if/when it's ever enabled
        factors.push_back({frac, {factors.back().second[0], 0., 0.}});
    }
    return factors;
}

double FractionalIntegralModel::get_integral_remainder(const Channel &channel) const {
    return (1. - this->at(channel)->get_value()) * _get_integral_remainder(_parent.get(), *_model, channel);
}

ProperFractionParameterD &FractionalIntegralModel::get_parameter_frac(const Channel &channel) const {
    return *(_map.at(channel));
}

ParamRefs &FractionalIntegralModel::get_parameters(ParamRefs &params, ParamFilter *filter) const {
    _model->get_parameters(params, filter);
    // Don't return the n_channels fixed frac=1 parameters at the end
    const size_t n_p_max = _data.size() - (this->is_final() ? this->get_channels().size() : 0);
    size_t i = 0;
    for (auto &p : _data) {
        if (i++ == n_p_max) break;
        insert_param_channel(p.first, *p.second, params, filter);
    }
    return params;
}

ParamCRefs &FractionalIntegralModel::get_parameters_const(ParamCRefs &params, ParamFilter *filter) const {
    _model->get_parameters_const(params, filter);
    const size_t n_p_max = _data.size() - (this->is_final() ? this->get_channels().size() : 0);
    size_t i = 0;
    for (const auto &p : _data) {
        if (i++ == n_p_max) break;
        insert_param_channel(p.first, *p.second, params, filter);
    }
    return params;
}

bool FractionalIntegralModel::is_final() const { return _is_final; }

size_t FractionalIntegralModel::size() const { return _data.size(); }

std::string FractionalIntegralModel::repr(bool name_keywords, std::string_view namespace_separator) const {
    std::string s = type_name_str<FractionalIntegralModel>(false, namespace_separator) + "("
                    + (name_keywords ? "data={" : "{");
    for (const auto &datum : _data) {
        s += datum.first.get().repr(name_keywords, namespace_separator) + ": "
             + datum.second->repr(name_keywords, namespace_separator) + ",";
    }
    return s + "})";
}

std::string FractionalIntegralModel::str() const {
    std::string s = type_name_str<FractionalIntegralModel>(true) + "(data={";
    for (const auto &datum : _data) {
        s += datum.first.get().str() + ": " + datum.second->str() + ",";
    }
    s += "}, model=" + _model->str() + ", is_final=" + std::to_string(_is_final) + ")";
    return s;
}

// Return a pointer to a registered FractionalIntegralModel, if it is one
// TODO: Verify if this is safe, or if dynamic_cast is better
// TODO: Consider a second constructor (but it would likely complicate pybind11'ing)
std::shared_ptr<const FractionalIntegralModel> FractionalIntegralModel::_find_parent(
        std::shared_ptr<const IntegralModel> model) {
    if (model == nullptr) return nullptr;
    auto parent = std::dynamic_pointer_cast<const FractionalIntegralModel>(model);
    if (parent == nullptr) return nullptr;
    bool found = FractionalIntegralModel::_registry.find(*parent) != _registry.end();
    return found ? parent : nullptr;
}

std::shared_ptr<FractionalIntegralModel> FractionalIntegralModel::make(
        std::optional<const Data> data, std::shared_ptr<const IntegralModel> model, bool is_final) {
    std::shared_ptr<FractionalIntegralModel> ptr
            = std::make_shared<Shared_enabler>(data, std::move(model), is_final);
    ;
    _registry.insert({*ptr, *model});
    FractionalIntegralModel::_registry_rev.insert({*model, ptr});
    return ptr;
}

const std::shared_ptr<const FractionalIntegralModel> FractionalIntegralModel::make_const(
        std::optional<const Data> data, std::shared_ptr<const IntegralModel> model, bool is_final) {
    return make(data, model, is_final);
}

}  // namespace lsst::gauss2d::fit
