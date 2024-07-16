#ifndef LSST_GAUSS2D_FIT_DATA_H
#define LSST_GAUSS2D_FIT_DATA_H

#include <memory>
#include <stdexcept>
#include <vector>

#include "channel.h"
#include "chromatic.h"
#include "observation.h"
#include "param_defs.h"
#include "parametric.h"

namespace lsst::gauss2d::fit {

/*
    A Data is a collection of Observations that a consistent Model can be
    generated for.
*/
/**
 * @brief A list of Observation instances that can be modelled.
 *
 * A Data is a list of Observation instances that can have associated Model
 * instances. Multiple Observation instances may have the same Channel, but one
 * should not include the same Observation multiple time.
 *
 * @tparam T The type of the Observation Image (usually float or double)
 * @tparam I The type of the Observation indices (usually size_t)
 * @tparam M The type of the Observation Mask (usually bool)
 */
template <typename T, typename I, typename M>
class Data : public Chromatic, public Parametric {
public:
    using Observation = lsst::gauss2d::fit::Observation<T, I, M>;
    using ObservationCRef = std::reference_wrapper<const Observation>;

    /**
     * Construct a Data instance
     *
     * @param observations The Observation pointers to include. Must not be null.
     */
    explicit Data(std::vector<std::shared_ptr<const Observation>> observations) {
        _observations.reserve(observations.size());
        _observation_ptrs.reserve(observations.size());

        for (const auto& observation : observations) {
            if (observation == nullptr) throw std::invalid_argument("Can't store null Observation");
            const auto& channel = observation->get_channel();
            if (_channels.find(channel) == _channels.end()) {
                _channels_ordered.push_back(channel);
                _channels.insert(channel);
            }
            _observations.push_back(ObservationCRef(*observation));
            _observation_ptrs.push_back(observation);
        }
    }

    inline auto at(size_t i) const { return _observations.at(i); }
    inline auto begin() const { return _observations.begin(); }
    inline auto end() const { return _observations.end(); }

    inline auto cbegin() const { return _observations.begin(); }
    inline auto cend() const { return _observations.end(); }

    std::vector<std::reference_wrapper<const Channel>> get_channels() const override {
        return _channels_ordered;
    }

    ParamRefs& get_parameters(ParamRefs& params, ParamFilter* filter = nullptr) const override {
        for (const Observation& exp : *this) exp.get_parameters(params, filter);
        return params;
    }
    ParamCRefs& get_parameters_const(ParamCRefs& params, ParamFilter* filter = nullptr) const override {
        for (auto exp_it = this->cbegin(); exp_it != this->cend(); ++exp_it) {
            (*exp_it).get().get_parameters_const(params, filter);
        }
        return params;
    }

    /// Get the number of member Observation
    size_t size() const { return _observations.size(); }

    std::string repr(bool name_keywords = false,
                     std::string_view namespace_separator = Object::CC_NAMESPACE_SEPARATOR) const override {
        std::string str = std::string("Data(") + (name_keywords ? "observations=[" : "[");
        for (auto exp_it = this->cbegin(); exp_it != this->cend(); ++exp_it) {
            str += (*exp_it).get().repr(name_keywords, namespace_separator) + ",";
        }
        str += "]);";
        return str;
    }
    std::string str() const override {
        std::string str = "Data(observations=[";
        for (auto exp_it = this->cbegin(); exp_it != this->cend(); ++exp_it) {
            str += (*exp_it).get().str() + ",";
        }
        str += "])";
        return str;
    }

private:
    // This could be unordered, but std::hash<std::string> won't take const strings
    std::set<std::reference_wrapper<const Channel>> _channels = {};
    std::vector<std::reference_wrapper<const Channel>> _channels_ordered = {};
    std::vector<std::shared_ptr<const Observation>> _observation_ptrs = {};
    std::vector<ObservationCRef> _observations = {};
};

}  // namespace lsst::gauss2d::fit

#endif
