#ifndef GAUSS2D_FIT_FRACTIONALINTEGRALMODEL_H
#define GAUSS2D_FIT_FRACTIONALINTEGRALMODEL_H

#include <map>
#include <memory>
#include <optional>
#include <vector>

#include "component.h"
#include "param_defs.h"
#include "param_filter.h"
#include "parameters.h"
#include "integralmodel.h"

namespace gauss2d
{
namespace fit
{

class FractionalGaussianIntegral;

/*
    FractionalIntegralModel is a model that returns a parameterized fraction
    of the flux of another model.

    FractionalIntegralModels can be chained, but each IntegralModel can only
    have one dependent FractionalIntegralModel. The registry enforces this.

    FractionalIntegralModel chains should end in a FractionalIntegralModel
    with a ProperFractionParameter fixed at 1. This is not yet enforced.

    One consideration for fitting is that for lengthy chains, a single
    ProperFractionParameter valued at 1 midway through the chain will set
    the integrals for the rest of the chain to zero, regardless of the
    values of their ProperFractionParameters. One should consider setting
    more restrictive upper limits on ProperFractionParameter values in 
    such longer chains.
*/
class FractionalIntegralModel : public IntegralModel {
public:
    typedef std::map<
        std::reference_wrapper<const Channel>,
        std::shared_ptr<ProperFractionParameter>,
        std::less<const Channel>
    > Data;

private:
    Data _data = {};

    // TODO: See if all raw pointers can be changed to reference_wrappers or weak_ptrs
    std::shared_ptr<const FractionalIntegralModel> _find_parent(std::shared_ptr<const IntegralModel> model);
    const std::shared_ptr<const IntegralModel> _model;
    std::shared_ptr<const FractionalIntegralModel> _parent;
    static inline std::map<
        std::reference_wrapper<const FractionalIntegralModel>, 
        std::reference_wrapper<const IntegralModel>
    > _registry = {};
    static inline std::map<
        std::reference_wrapper<const IntegralModel>,
        std::weak_ptr<FractionalIntegralModel>
    > _registry_rev = {};

    struct Shared_enabler;

    bool _is_final;

    // not giving a nullptr default data_in because the map needs to match the model's channels
    FractionalIntegralModel(std::optional<const Data> data, std::shared_ptr<const IntegralModel> model,
        bool is_final);

public:
    std::shared_ptr<ProperFractionParameter> at(const Channel & channel);
    std::shared_ptr<const ProperFractionParameter> at(const Channel & channel) const;

    typename Data::iterator begin() noexcept;
    typename Data::const_iterator cbegin() const noexcept;

    typename Data::iterator end() noexcept;
    typename Data::const_iterator cend() const noexcept;

    // Find the FractionalIntegralModel that depends on a given IntegralModel, if any
    static std::shared_ptr<FractionalIntegralModel> find_model(const IntegralModel & model) {
        const auto found = _registry_rev.find(model);
        return (found == _registry_rev.end()) ? nullptr : (*found).second.lock();
    }

    std::set<std::reference_wrapper<const Channel>> get_channels() const override;
    const IntegralModel & get_parent_model() const;
    double get_integral(const Channel & channel) const override;
        std::vector<std::pair<ParamBaseCRef, extra_param_factor_values>> get_integral_derivative_factors(
        const Channel & channel) const override;
    double get_integral_remainder(const Channel & channel) const;

    ProperFractionParameter & get_parameter_frac(const Channel & channel) const;
    ParamRefs & get_parameters(ParamRefs & params, ParamFilter * filter = nullptr) const override;
    ParamCRefs & get_parameters_const(ParamCRefs & params, ParamFilter * filter = nullptr) const override;
    
    bool is_final() const;

    static std::shared_ptr<FractionalIntegralModel> make(
        std::optional<const Data> data,
        std::shared_ptr<const IntegralModel> model,
        bool is_final=false
    );
    static const std::shared_ptr<const FractionalIntegralModel> make_const(
        std::optional<const Data> data,
        std::shared_ptr<const IntegralModel> model,
        bool is_final=false
    );

    size_t size() const;

    std::string str() const override;

    FractionalIntegralModel (const FractionalIntegralModel&) = delete;
    FractionalIntegralModel& operator= (const FractionalIntegralModel&) = delete;

    ~FractionalIntegralModel();
};

} // namespace fit
} // namespace gauss2d

#endif
