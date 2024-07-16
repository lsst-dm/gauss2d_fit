#ifdef LSST_GAUSS2D_FIT_HAS_GSL

#include <stdexcept>
#include <vector>

#include "lsst/gauss2d/type_name.h"

#include "lsst/gauss2d/fit/gslinterpolator.h"

namespace lsst::gauss2d::fit {

void _check_idx(size_t idx, size_t size) {
    if (!(idx < size)) {
        throw std::out_of_range("idx=" + std::to_string(idx) + " !< size=" + std::to_string(size));
    }
}

GSLInterpolator::GSLInterpolator(std::vector<double> x, std::vector<double> y, const InterpType interp_type)
        : _interp_type(interp_type), _n_knots(x.size()), _x(x), _y(y) {
    if (!(_n_knots > 0)) {
        throw std::invalid_argument("x.size()=" + std::to_string(_n_knots) + " !> 0");
    }
    if (_y.size() != _n_knots) {
        throw std::invalid_argument("y.size()=" + std::to_string(y.size())
                                    + " != x.size()=" + std::to_string(_n_knots));
    }
    _acc = gsl_interp_accel_alloc();
    _spline = gsl_spline_alloc(GSLInterpTypes.at(interp_type), _n_knots);
    gsl_spline_init(_spline, &x[0], &y[0], _n_knots);
}

GSLInterpolator::~GSLInterpolator() {
    gsl_spline_free(_spline);
    gsl_interp_accel_free(_acc);
}

const InterpType GSLInterpolator::get_interp_type() const { return _interp_type; }

double GSLInterpolator::get_knot_x(size_t idx) const {
    _check_idx(idx, this->_n_knots);
    return _x[idx];
}
double GSLInterpolator::get_knot_y(size_t idx) const {
    _check_idx(idx, this->_n_knots);
    return _y[idx];
}

double GSLInterpolator::eval(double x) const { return gsl_spline_eval(_spline, x, _acc); }

double GSLInterpolator::eval_deriv(double x) const { return gsl_spline_eval_deriv(_spline, x, _acc); }

size_t GSLInterpolator::size() const { return _n_knots; }

std::string GSLInterpolator::repr(bool name_keywords, std::string_view namespace_separator) const {
    return std::string("GSLInterpolator(") + (name_keywords ? "n_knots=" : "") + std::to_string(_n_knots)
           + ")";
}

std::string GSLInterpolator::str() const {
    return "GSLInterpolator(n_knots=" + std::to_string(_n_knots) + ")";
}
}  // namespace lsst::gauss2d::fit

#endif  // #ifdef LSST_GAUSS2D_FIT_HAS_GSL
