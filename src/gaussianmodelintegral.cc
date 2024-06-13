#include "lsst/gauss2d/type_name.h"

#include "lsst/gauss2d/fit/channel.h"
#include "lsst/gauss2d/fit/gaussianmodelintegral.h"
#include "lsst/gauss2d/fit/integralmodel.h"

namespace lsst::gauss2d::fit {

GaussianModelIntegral::GaussianModelIntegral(const Channel& channel,
                                             const std::shared_ptr<const IntegralModel> integralmodel)
        : _channel(channel), _integralmodel(std::move(integralmodel)) {
    if (_integralmodel == nullptr)
        throw std::invalid_argument("GaussianModelIntegral integralmodel can't be null");
}
GaussianModelIntegral::~GaussianModelIntegral(){};

double GaussianModelIntegral::get_value() const { return _integralmodel->get_integral(_channel); }
void GaussianModelIntegral::set_value(double value) {
    throw std::runtime_error("Can't set_value on GaussianModelIntegral");
}

std::string GaussianModelIntegral::repr(bool name_keywords, std::string_view namespace_separator) const {
    return type_name_str<GaussianModelIntegral>(false, namespace_separator) + "("
           + (name_keywords ? "channel=" : "") + _channel.repr(name_keywords, namespace_separator) + ", "
           + (name_keywords ? "integralmodel=" : "")
           + _integralmodel->repr(name_keywords, namespace_separator) + ")";
}

std::string GaussianModelIntegral::str() const {
    return type_name_str<GaussianModelIntegral>(true) + "(channel=" + _channel.str()
           + ", integralmodel=" + _integralmodel->str() + ")";
}

}  // namespace lsst::gauss2d::fit
