#ifndef GAUSS2D_FIT_LINEARINTEGRALMODEL_H
#define GAUSS2D_FIT_LINEARINTEGRALMODEL_H

#include <map>
#include <memory>

#include "channel.h"
#include "parameters.h"
#include "integralmodel.h"

namespace gauss2d
{
namespace fit
{

class LinearIntegralModel : public IntegralModel
{
public:
    typedef std::map<
        std::reference_wrapper<const Channel>,
        std::shared_ptr<IntegralParameter>,
        std::less<const Channel>
    > Data;

private:
    Data _data = {};
    struct Shared_enabler;

public:
    std::set<std::reference_wrapper<const Channel>> get_channels() const override;
    double get_integral(const Channel & channel) const override;
    
    std::shared_ptr<IntegralParameter> at(const Channel & channel);
    std::shared_ptr<const IntegralParameter> at(const Channel & channel) const;

    typename Data::iterator begin() noexcept;
    typename Data::const_iterator cbegin() const noexcept;

    typename Data::iterator end() noexcept;
    typename Data::const_iterator cend() const noexcept;

    ParamRefs & get_parameters(ParamRefs & params, ParamFilter * filter = nullptr) const override;
    ParamCRefs & get_parameters_const(ParamCRefs & params, ParamFilter * filter = nullptr) const override;

    size_t size() const;

    std::string str() const override;

    const bool operator < ( const IntegralModel &m ) const { return &(*this) < &m; };
    const bool operator == ( const IntegralModel &m ) const { return &(*this) == &m; };
    const bool operator != ( const IntegralModel &m ) const { return &(*this) != &m; };


    // not giving a nullptr default because users should explicitly use the null Channel
    LinearIntegralModel(const Data * data_in);
    ~LinearIntegralModel();
};


} // namespace fit
} // namespace gauss2d

#endif