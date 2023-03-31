#ifndef GAUSS2D_FIT_DATA_H
#define GAUSS2D_FIT_DATA_H

#include <memory>
#include <stdexcept>
#include <vector>

#include "observation.h"
#include "param_defs.h"
#include "parametric.h"

namespace gauss2d::fit
{

/*
    A Data is a collection of Observations that a consistent Model can be
    generated for.
*/
template <typename T, typename I, typename M>
class Data : public Parametric
{
public:
    using Observation = gauss2d::fit::Observation<T, I, M>;
    using ObservationCRef = std::reference_wrapper<const Observation>;

private:
    std::set<std::reference_wrapper<const Channel>> _channels {};
    std::vector<std::shared_ptr<const Observation>> _observation_ptrs {};
    std::vector<ObservationCRef> _observations {};

public:
    inline auto at(size_t i) const { return _observations.at(i); }
    inline auto begin() const { return _observations.begin(); }
    inline auto end() const { return _observations.end(); }

    inline auto cbegin() const { return _observations.begin(); }
    inline auto cend() const { return _observations.end(); }

    std::set<std::reference_wrapper<const Channel>> get_channels() const {
        return _channels;
    }

    ParamRefs & get_parameters(ParamRefs & params, ParamFilter * filter = nullptr) const override {
        for(const Observation & exp: *this) exp.get_parameters(params, filter);
        return params;
    }
    ParamCRefs & get_parameters_const(ParamCRefs & params, ParamFilter * filter = nullptr) const override {
        for(auto exp_it = this->cbegin(); exp_it != this->cend(); ++exp_it) {
            (*exp_it).get().get_parameters_const(params, filter);
        }
        return params;
    }

    size_t size() const { return _observations.size(); }

    std::string repr(bool name_keywords = false) const override {
        std::string str = std::string("Data(") + (name_keywords ? "observations=[" : "[");
        for(auto exp_it = this->cbegin(); exp_it != this->cend(); ++exp_it) {
            str += (*exp_it).get().repr(name_keywords) + ",";
        }
        str += "]);";
        return str;
    }
    std::string str() const override {
        std::string str = "Data(observations=[";
        for(auto exp_it = this->cbegin(); exp_it != this->cend(); ++exp_it) {
            str += (*exp_it).get().str() + ",";
        }
        str += "]);";
        return str;
    }

    explicit Data(std::vector<std::shared_ptr<const Observation>> observations) {
        _observations.reserve(observations.size());
        _observation_ptrs.reserve(observations.size());
        for(const auto & observation: observations) {
            if(observation == nullptr) throw std::invalid_argument("Can't store null Observation");
            _channels.insert(observation->get_channel());
            _observations.push_back(ObservationCRef(*observation));
            _observation_ptrs.push_back(std::move(observation));
        }
    }
};

} // namespace gauss2d::fit

#endif
