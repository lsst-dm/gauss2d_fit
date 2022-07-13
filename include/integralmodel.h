#ifndef GAUSS2D_FIT_INTEGRALMODEL_H
#define GAUSS2D_FIT_INTEGRALMODEL_H

#include <memory>
#include <set>

#include "channel.h"
#include "parametric.h"

namespace gauss2d
{
namespace fit
{

/*
    An IntegralModel described the total integrated surface density of a
    2D distribution. It is intended to defined a Component but can
    be used otherwise.

    It is generally not intended for IntegralModels to be shared amonst
    Components or otherwise, but this is not guaranteed to be enforced.
*/
class IntegralModel : public Parametric
{
public:
    virtual std::set<std::reference_wrapper<const Channel>> get_channels() const = 0;
    virtual double get_integral(const Channel & channel) const = 0;
};

inline bool operator < ( const IntegralModel &lhs, const IntegralModel &rhs) { return &lhs < &rhs; }
//const bool operator == ( const IntegralModel &m ) const { return &(*this) == &m; };
//const bool operator != ( const IntegralModel &m ) const { return &(*this) != &m; };


} // namespace fit
} // namespace gauss2d

#endif