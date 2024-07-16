#ifndef LSST_GAUSS2D_FIT_FRACTIONALINTEGRALMODEL_H
#define LSST_GAUSS2D_FIT_FRACTIONALINTEGRALMODEL_H

#include <map>
#include <memory>
#include <optional>
#include <vector>

#include "component.h"
#include "param_defs.h"
#include "param_filter.h"
#include "parameters.h"
#include "integralmodel.h"

namespace lsst::gauss2d::fit {
/**
 * @brief An IntegralModel that returns a Parameter-dependent fraction
 * of the flux of another IntegralModel.
 *
 * FractionalIntegralModel instances can be chained but are registered
 * to enforce having at most one dependent FractionalIntegralModel. Chains
 * should end in a single fixed ProperFractionParameter of value 1; this
 * is not yet fully enforced.
 *
 * @note For longer chains, a single ProperFractionParameter of value 1
 * will render all downstream models zero-valued, regardless of their own
 * ProperFractionParameter values. Consider setting Limits to prevent such
 * situations if they become problematic.
 */
class FractionalIntegralModel : public IntegralModel {
public:
    typedef std::pair<std::reference_wrapper<const Channel>, std::shared_ptr<ProperFractionParameterD>>
            ChannelIntegralParameterD;
    typedef std::vector<ChannelIntegralParameterD> Data;

    explicit FractionalIntegralModel(const FractionalIntegralModel&) = delete;
    FractionalIntegralModel& operator=(const FractionalIntegralModel&) = delete;

    ~FractionalIntegralModel();

    std::shared_ptr<ProperFractionParameterD> at(const Channel& channel);
    std::shared_ptr<const ProperFractionParameterD> at(const Channel& channel) const;

    typename Data::iterator begin() noexcept;
    typename Data::const_iterator cbegin() const noexcept;

    typename Data::iterator end() noexcept;
    typename Data::const_iterator cend() const noexcept;

    /**
     * Find the FractionalIntegralModel that depends on a given IntegralModel, if any.
     *
     * @param model The IntegralModel to search for
     * @return The FractionalIntegralModel that depends on model, or nullptr if none
     */
    static std::shared_ptr<FractionalIntegralModel> find_model(const IntegralModel& model) {
        const auto found = _registry_rev.find(model);
        return (found == _registry_rev.end()) ? nullptr : (*found).second.lock();
    }

    std::vector<std::reference_wrapper<const Channel>> get_channels() const override;
    const IntegralModel& get_parent_model() const;
    double get_integral(const Channel& channel) const override;
    std::vector<std::pair<ParamBaseCRef, ExtraParamFactorValues>> get_integral_derivative_factors(
            const Channel& channel) const override;
    double get_integral_remainder(const Channel& channel) const;

    ProperFractionParameterD& get_parameter_frac(const Channel& channel) const;
    ParamRefs& get_parameters(ParamRefs& params, ParamFilter* filter = nullptr) const override;
    ParamCRefs& get_parameters_const(ParamCRefs& params, ParamFilter* filter = nullptr) const override;

    bool is_final() const;

    /**
     * Construct a FractionalIntegralModel and add to registry
     *
     * @param data The map of ProperFractionParameter instances for each channel.
     * @param model The IntegralModel that ProperFractionParameter instances multiply by.
     * @param is_final Whether these fractions are the last in the chain.
     *
     * @note Input parameters are validated and will throw exceptions if invalid, e.g.
     *       if data and model do not cover the same channels.
     *
     * @return A new FractionalIntegralModel instance
     */
    static std::shared_ptr<FractionalIntegralModel> make(std::optional<const Data> data,
                                                         std::shared_ptr<const IntegralModel> model,
                                                         bool is_final = false);
    static const std::shared_ptr<const FractionalIntegralModel> make_const(
            std::optional<const Data> data, std::shared_ptr<const IntegralModel> model,
            bool is_final = false);

    size_t size() const;

    std::string repr(bool name_keywords = false,
                     std::string_view namespace_separator = Object::CC_NAMESPACE_SEPARATOR) const override;
    std::string str() const override;

private:
    Data _data = {};
    // This could be unordered, but std::hash<std::string> won't take const strings
    // (see also linearintegralmodel.h)
    std::map<std::reference_wrapper<const Channel>, std::shared_ptr<ProperFractionParameterD>> _map = {};

    // TODO: See if all raw pointers can be changed to reference_wrappers or weak_ptrs
    std::shared_ptr<const FractionalIntegralModel> _find_parent(std::shared_ptr<const IntegralModel> model);
    const std::shared_ptr<const IntegralModel> _model;
    std::shared_ptr<const FractionalIntegralModel> _parent;
    static inline std::map<std::reference_wrapper<const FractionalIntegralModel>,
                           std::reference_wrapper<const IntegralModel>>
            _registry = {};
    static inline std::map<std::reference_wrapper<const IntegralModel>,
                           std::weak_ptr<FractionalIntegralModel>>
            _registry_rev = {};

    struct Shared_enabler;

    bool _is_final;

    // not giving a nullptr default data_in because the map needs to match the model's channels
    FractionalIntegralModel(std::optional<const Data> data, std::shared_ptr<const IntegralModel> model,
                            bool is_final);
};

}  // namespace lsst::gauss2d::fit

#endif
