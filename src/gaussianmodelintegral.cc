#include "gaussianmodelintegral.h"

#include "channel.h"
#include "integralmodel.h"

namespace gauss2d
{
namespace fit
{

double GaussianModelIntegral::get_value() const {
    return _integralmodel->get_integral(_channel);
}
void GaussianModelIntegral::set_value(double value) {
    throw std::runtime_error("Can't set_value on GaussianModelIntegral");
}

std::string GaussianModelIntegral::repr(bool name_keywords) const {
    return std::string("GaussianModelIntegral(")
        + (name_keywords ? "channel=" : "") + _channel.repr(name_keywords) + ", "
        + (name_keywords ? "integralmodel=" : "") + _integralmodel->repr(name_keywords) + ")";
}

std::string GaussianModelIntegral::str() const {
    return "GaussianModelIntegral(channel=" + _channel.str()
        + ", integralmodel=" + _integralmodel->str() + ")";
}

GaussianModelIntegral::GaussianModelIntegral(
    const Channel & channel, const std::shared_ptr<const IntegralModel> integralmodel
) : _channel(channel), _integralmodel(std::move(integralmodel))
{
    if(_integralmodel == nullptr) throw std::invalid_argument("GaussianModelIntegral integralmodel can't be null");
}
GaussianModelIntegral::~GaussianModelIntegral() {};

} // namespace fit
} // namespace gauss2d
