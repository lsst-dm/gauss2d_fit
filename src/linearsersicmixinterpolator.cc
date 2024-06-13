#include <algorithm>
#include <stdexcept>
#include <string>
#include <vector>

#include "lsst/gauss2d/to_string.h"
#include "lsst/gauss2d/type_name.h"

#include "lsst/gauss2d/fit/linearsersicmixinterpolator.h"

namespace lsst::gauss2d::fit {

LinearSersicMixInterpolator::LinearSersicMixInterpolator(unsigned short order)
        : _order(order),
          _knots(get_sersic_mix_knots(order)),
          _sersicindex_min(_knots[0].sersicindex),
          _sersicindex_max(_knots.back().sersicindex) {}

LinearSersicMixInterpolator::~LinearSersicMixInterpolator() {}

std::vector<IntegralSize> LinearSersicMixInterpolator::get_integralsizes(double sersicindex) const {
    if (!((sersicindex >= _sersicindex_min) && (sersicindex <= _sersicindex_max))) {
        throw std::invalid_argument("sersicindex=" + to_string_float(sersicindex)
                                    + " !(>=" + to_string_float(_sersicindex_min)
                                    + "&& <=" + to_string_float(_sersicindex_max));
    }

    if (sersicindex == _sersicindex_min)
        return _knots[0].values;
    else if (sersicindex == _sersicindex_max)
        return _knots.back().values;

    auto found = std::lower_bound(_knots.begin(), _knots.end(), sersicindex);
    auto high = *found;
    auto low = *(--found);

    if (sersicindex == high.sersicindex) return high.values;
    double frac_low = (high.sersicindex - sersicindex) / (high.sersicindex - low.sersicindex);
    if (!((frac_low >= 0) && (frac_low <= 1))) {
        throw std::logic_error("Got invalid frac_low=" + std::to_string(frac_low)
                               + " with n, lo, hi=" + to_string_float(sersicindex) + ","
                               + to_string_float(low.sersicindex) + "," + to_string_float(high.sersicindex));
    }
    double frac_high = 1 - frac_low;

    std::vector<IntegralSize> result;
    result.reserve(_order);
    for (size_t i = 0; i < _order; ++i) {
        result.push_back(IntegralSize(frac_low * low.values[i].integral + frac_high * high.values[i].integral,
                                      frac_low * low.values[i].sigma + frac_high * high.values[i].sigma));
    }

    return result;
}

std::vector<IntegralSize> LinearSersicMixInterpolator::get_integralsizes_derivs(double sersicindex) const {
    if (!((sersicindex >= _sersicindex_min) && (sersicindex <= _sersicindex_max))) {
        throw std::invalid_argument("sersicindex=" + to_string_float(sersicindex)
                                    + " !(>=" + to_string_float(_sersicindex_min)
                                    + "&& <=" + to_string_float(_sersicindex_max));
    }

    auto found = sersicindex == _sersicindex_min
                         ? ++_knots.begin()
                         : ((sersicindex == _sersicindex_max)
                                    ? --_knots.end()
                                    : std::upper_bound(_knots.begin(), _knots.end(), sersicindex));
    auto& high = *found;
    auto& low = *(--found);

    double dn_inv = 1. / (high.sersicindex - low.sersicindex);

    std::vector<IntegralSize> result;
    result.reserve(_order);
    for (size_t i = 0; i < _order; ++i) {
        result.push_back(IntegralSize((high.values[i].integral - low.values[i].integral) * dn_inv,
                                      (high.values[i].sigma - low.values[i].sigma) * dn_inv));
    }

    return result;
}

InterpType LinearSersicMixInterpolator::get_interptype() const { return InterpType::linear; }

const std::vector<SersicMixValues>& LinearSersicMixInterpolator::get_knots() const { return _knots; }

unsigned short LinearSersicMixInterpolator::get_order() const { return _order; }

double LinearSersicMixInterpolator::get_sersicindex_min() const { return _sersicindex_min; }

double LinearSersicMixInterpolator::get_sersicindex_max() const { return _sersicindex_max; }

std::string LinearSersicMixInterpolator::repr(bool name_keywords,
                                              std::string_view namespace_separator) const {
    return type_name_str<LinearSersicMixInterpolator>(false, namespace_separator) + "("
           + (name_keywords ? "order=" : "") + std::to_string(_order) + ")";
}

std::string LinearSersicMixInterpolator::str() const {
    return type_name_str<LinearSersicMixInterpolator>(true) + "(order=" + std::to_string(_order) + ")";
}

}  // namespace lsst::gauss2d::fit
