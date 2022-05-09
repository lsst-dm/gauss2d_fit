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

    // not giving a nullptr default data_in because the map needs to match the model's channels
    FractionalIntegralModel(std::optional<const Data> data, std::shared_ptr<const IntegralModel> model);

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
    double get_integral_remainder(const Channel & channel) const;

    ParamRefs & get_parameters(ParamRefs & params, ParamFilter * filter = nullptr) const override;
    ParamCRefs & get_parameters_const(ParamCRefs & params, ParamFilter * filter = nullptr) const override;
    
    static std::shared_ptr<FractionalIntegralModel> make(
        std::optional<const Data> data,
        std::shared_ptr<const IntegralModel> model
    );
    static const std::shared_ptr<const FractionalIntegralModel> make_const(
        std::optional<const Data> data,
        std::shared_ptr<const IntegralModel> model
    );

    size_t size() const;

    std::string str() const override;

    FractionalIntegralModel (const FractionalIntegralModel&) = delete;
    FractionalIntegralModel& operator= (const FractionalIntegralModel&) = delete;

    ~FractionalIntegralModel();
};

class FractionalGaussianIntegral : public GaussianIntegral, public Parametric
{
private:
    const std::shared_ptr<Channel> _channel;
    std::shared_ptr<FractionalIntegralModel> _model;

public:
    ParamRefs & get_parameters(ParamRefs & params, ParamFilter * filter = nullptr) const override;
    ParamCRefs & get_parameters_const(ParamCRefs & params, ParamFilter * filter = nullptr) const override;

    double get_value() const override;
    void set_value(double value) override;

    std::string str() const override;

    FractionalGaussianIntegral(
        const std::shared_ptr<Channel> _channel,
        std::shared_ptr<FractionalIntegralModel> param_frac
    );
    ~FractionalGaussianIntegral() {};
};


} // namespace fit
} // namespace gauss2d

#endif