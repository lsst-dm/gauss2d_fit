#ifndef GAUSS2D_FIT_LINEARINTEGRALMODEL_H
#define GAUSS2D_FIT_LINEARINTEGRALMODEL_H

#include <map>
#include <memory>

#include "channel.h"
#include "parameters.h"
#include "integralmodel.h"

namespace gauss2d::fit {

/**
 * An IntegralModel with a single linear IntegralParameter per Channel.
 */
class LinearIntegralModel : public IntegralModel {
public:
    typedef std::pair<std::reference_wrapper<const Channel>, std::shared_ptr<IntegralParameter>>
            ChannelIntegralParameter;
    typedef std::vector<ChannelIntegralParameter> Data;

private:
    Data _data = {};
    // This could be unordered, but std::hash<std::string> won't take const strings
    // ... and it doesn't seem to be worth the effort to work around
    std::map<std::reference_wrapper<const Channel>, std::shared_ptr<IntegralParameter>> _map = {};
    struct Shared_enabler;

public:
    /// Get the IntegralParameter for the given Channel
    std::shared_ptr<IntegralParameter> at(const Channel &channel);
    /// Get the (const) IntegralParameter for the given Channel
    std::shared_ptr<const IntegralParameter> at(const Channel &channel) const;

    typename Data::iterator begin() noexcept;
    typename Data::const_iterator cbegin() const noexcept;

    typename Data::iterator end() noexcept;
    typename Data::const_iterator cend() const noexcept;

    std::vector<std::reference_wrapper<const Channel>> get_channels() const override;
    double get_integral(const Channel &channel) const override;
    std::vector<std::pair<ParamBaseCRef, ExtraParamFactorValues>> get_integral_derivative_factors(
            const Channel &channel) const override;

    ParamRefs &get_parameters(ParamRefs &params, ParamFilter *filter = nullptr) const override;
    ParamCRefs &get_parameters_const(ParamCRefs &params, ParamFilter *filter = nullptr) const override;

    /// Return the size of Data (number of Channel/IntegralParameter instances)
    size_t size() const;

    std::string repr(bool name_keywords = false) const override;
    std::string str() const override;

    const bool operator<(const IntegralModel &m) const { return &(*this) < &m; };
    const bool operator==(const IntegralModel &m) const { return &(*this) == &m; };
    const bool operator!=(const IntegralModel &m) const { return &(*this) != &m; };

    /**
     * Construct a LinearIntegralModel from input Data.
     *
     * @param data_in A map of IntegralParameter shared_ptr to move for each Channel.
     *
     * @note No default initialization is provided, so data_in must not be null.
     */
    LinearIntegralModel(const Data *data_in);
    ~LinearIntegralModel();
};

}  // namespace gauss2d::fit

#endif
