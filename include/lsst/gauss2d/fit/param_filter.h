#ifndef LSST_GAUSS2D_FIT_PARAM_FILTER_H
#define LSST_GAUSS2D_FIT_PARAM_FILTER_H

#include <cstddef>
#include <optional>

#include "channel.h"
#include "param_defs.h"

namespace g2f = lsst::gauss2d::fit;

namespace lsst::gauss2d::fit{

/**
 * Options for filtering Parameter instances. Defaults select all.
 */
struct ParamFilter {
    bool fixed = true;
    bool free = true;
    bool linear = true;
    bool nonlinear = true;
    std::optional<std::reference_wrapper<const Channel>> channel = std::nullopt;
};

/**
 * Add a Parameter to a vector thereof, if it meets conditions.
 *
 * @tparam t The type of the Parameter class (usually double or float).
 * @param param The Parameter to add (or not).
 * @param params The vector to add to.
 * @param filter The filter to apply.
 */
template <typename t>
inline void insert_param(g2f::ParamBase& param, t& params, ParamFilter* filter = nullptr) {
    bool is_fixed = param.get_fixed();
    bool is_linear = param.get_linear();
    if ((filter == nullptr)
        || (((is_fixed == filter->fixed) || (is_fixed != filter->free))
            && ((is_linear == filter->linear) || (is_linear != filter->nonlinear)))) {
        params.insert(params.end(), param);
    }
}

/// Same as insert_param, but only inserting if filter has a null or matching Channel.
template <typename t>
void insert_param_channel(const lsst::gauss2d::fit::Channel& channel, g2f::ParamBase& param, t& params,
                          ParamFilter* filter) {
    if ((filter == nullptr) || (filter->channel == std::nullopt) || (filter->channel == channel)) {
        insert_param(param, params, filter);
    }
}

/// Same as insert_param, but inserting from a (ref to a) collection of Parameter instances.
template <typename t>
inline void insert_params_ref(const t& params_in, t& params_out, ParamFilter* filter = nullptr) {
    if (filter != nullptr)
        for (auto& p : params_in) insert_param<t>(p.get(), params_out, filter);
    else
        for (auto& p : params_in) params_out.push_back(p);
}

/// Same as insert_param, but inserting from a collection of Parameter instances.
template <typename t>
inline void insert_params(const t params_in, t& params_out, ParamFilter* filter = nullptr) {
    insert_params_ref<t>(params_in, params_out, filter);
}

}  // namespace lsst::gauss2d::fit
#endif  // GAUSS2D_FIT_PARAM_FILTER_H
