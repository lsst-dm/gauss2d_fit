#ifdef LSST_GAUSS2D_FIT_HAS_GSL

#include <algorithm>
#include <stdexcept>
#include <string>
#include <vector>

#include "lsst/gauss2d/to_string.h"
#include "lsst/gauss2d/type_name.h"

#include "lsst/gauss2d/fit/gslsersicmixinterpolator.h"
#include "lsst/gauss2d/fit/util.h"

namespace lsst::gauss2d::fit {

GSLSersicMixInterpolator::GSLSersicMixInterpolator(unsigned short order, const InterpType interp_type_)
        : _interp_type(interp_type_),
          _order(order),
          _knots(get_sersic_mix_knots(order)),
          _sersicindex_min(_knots[0].sersicindex),
          _sersicindex_max(_knots.back().sersicindex) {
    _interps.reserve(order);
    std::vector<double> sersics;
    std::vector<std::vector<double>> integrals;
    std::vector<std::vector<double>> sigmas;
    const size_t n_knots = _knots.size();
    sersics.reserve(n_knots);
    integrals.resize(order);
    sigmas.resize(order);
    for (size_t iord = 0; iord < _order; ++iord) {
        integrals[iord].reserve(n_knots);
        sigmas[iord].reserve(n_knots);
    }
    for (const auto& knot : _knots) {
        sersics.push_back(knot.sersicindex);
        for (size_t iord = 0; iord < _order; ++iord) {
            const auto& integralsize = knot.values[iord];
            integrals[iord].push_back(integralsize.integral);
            sigmas[iord].push_back(integralsize.sigma);
        }
    }
    for (size_t iord = 0; iord < _order; ++iord) {
        _interps.push_back({std::make_unique<GSLInterpolator>(sersics, integrals[iord], _interp_type),
                            std::make_unique<GSLInterpolator>(sersics, sigmas[iord], _interp_type)});
    }
}

GSLSersicMixInterpolator::~GSLSersicMixInterpolator() {}

double GSLSersicMixInterpolator::get_final_correction() const { return _final_correction; }

std::vector<IntegralSize> GSLSersicMixInterpolator::get_integralsizes(double sersicindex) const {
    if (!((sersicindex >= _sersicindex_min) && (sersicindex <= _sersicindex_max))) {
        throw std::invalid_argument("sersicindex=" + to_string_float(sersicindex)
                                    + " not between min=" + to_string_float(_sersicindex_min)
                                    + " and max=" + to_string_float(_sersicindex_max));
    }

    std::vector<IntegralSize> result;
    result.reserve(_order);
    double sum = 0;
    size_t max_ord = _order - correct_final_integral;
    for (size_t i = 0; i < max_ord; ++i) {
        const auto& interp = _interps[i];
        double integral = interp.first->eval(sersicindex);
        sum += integral;
        result.push_back(IntegralSize(integral, interp.second->eval(sersicindex)));
    }
    if (correct_final_integral) {
        const auto& interp = _interps[_order - 1];
        double integral = 1.0 - sum;
        double integral_interp = interp.first->eval(sersicindex);
        bool interp_eq_zero = integral_interp == 0;
        _final_correction = integral * (!interp_eq_zero) / (integral_interp + interp_eq_zero);
        if (!(integral >= 0)) {
            throw std::logic_error(this->str() + ".get_integralsizes(" + to_string_float(sersicindex)
                                   + ") got correction=" + to_string_float(_final_correction)
                                   + " from integral=" + to_string_float(integral) + " !>=0)");
        }
        result.push_back(IntegralSize(integral, interp.second->eval(sersicindex)));
    } else {
        // An overly cautious reset of the value since the option is mutable.
        _final_correction = 1.;
    }

    return result;
}

std::vector<IntegralSize> GSLSersicMixInterpolator::get_integralsizes_derivs(double sersicindex) const {
    if (!((sersicindex >= _sersicindex_min) && (sersicindex <= _sersicindex_max))) {
        throw std::invalid_argument("sersicindex=" + to_string_float(sersicindex)
                                    + " not between min=" + to_string_float(_sersicindex_min)
                                    + " and max=" + to_string_float(_sersicindex_max));
    }

    std::vector<IntegralSize> result;
    result.reserve(_order);
    size_t max_ord = _order - correct_final_integral;
    for (size_t i = 0; i < max_ord; ++i) {
        const auto& interp = _interps[i];
        result.push_back(
                IntegralSize(interp.first->eval_deriv(sersicindex), interp.second->eval_deriv(sersicindex)));
    }
    if (correct_final_integral) {
        const auto& interp = _interps[max_ord];
        result.push_back(IntegralSize(correct_final_integral * interp.first->eval_deriv(sersicindex),
                                      interp.second->eval_deriv(sersicindex)));
    } else {
        _final_correction = 1.;
    }

    return result;
}

InterpType GSLSersicMixInterpolator::get_interptype() const { return _interp_type; }

const std::vector<SersicMixValues>& GSLSersicMixInterpolator::get_knots() const { return _knots; }

unsigned short GSLSersicMixInterpolator::get_order() const { return _order; }

double GSLSersicMixInterpolator::get_sersicindex_min() const { return _sersicindex_min; }
double GSLSersicMixInterpolator::get_sersicindex_max() const { return _sersicindex_max; }

std::string GSLSersicMixInterpolator::repr(bool name_keywords, std::string_view namespace_separator) const {
    return std::string("GSLSersicMixInterpolator(") + (name_keywords ? "order=" : "") + std::to_string(_order)
           + ")";
}

std::string GSLSersicMixInterpolator::str() const {
    return "GSLSersicMixInterpolator(order=" + std::to_string(_order) + ")";
}
}  // namespace lsst::gauss2d::fit

#endif  // #ifdef LSST_GAUSS2D_FIT_HAS_GSL
