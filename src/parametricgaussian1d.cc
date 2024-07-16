#include "lsst/gauss2d/type_name.h"

#include "lsst/gauss2d/fit/parameters.h"
#include "lsst/gauss2d/fit/parametricgaussian1d.h"

namespace lsst::gauss2d::fit {

ParametricGaussian1D::ParametricGaussian1D(std::shared_ptr<MeanParameterD> mean,
                                           std::shared_ptr<StdDevParameterD> stddev)
        : _mean(std::move(mean)), _stddev(std::move(stddev)) {
    if (_mean == nullptr) _mean = std::make_shared<MeanParameterD>();
    if (_stddev == nullptr) _stddev = std::make_shared<StdDevParameterD>();
}

ParametricGaussian1D::~ParametricGaussian1D() {}

double ParametricGaussian1D::get_mean() const { return this->_mean->get_value(); }

double ParametricGaussian1D::get_stddev() const { return this->_stddev->get_value(); }

MeanParameterD &ParametricGaussian1D::get_mean_parameter() const { return *this->_mean; }

StdDevParameterD &ParametricGaussian1D::get_stddev_parameter() const { return *this->_stddev; }

void ParametricGaussian1D::set_mean(double value) { this->_mean->set_value(value); }

void ParametricGaussian1D::set_stddev(double value) { this->_stddev->set_value(value); }

std::string ParametricGaussian1D::repr(bool name_keywords, std::string_view namespace_separator) const {
    return type_name_str<ParametricGaussian1D>(false, namespace_separator) + "("
           + (name_keywords ? "mean=" : "") + _mean->repr(name_keywords, namespace_separator)
           + (name_keywords ? ", stddev=" : "") + _stddev->repr(name_keywords, namespace_separator) + ")";
}

std::string ParametricGaussian1D::str() const {
    return type_name_str<ParametricGaussian1D>(true) + "(mean=" + _mean->str() + ", stddev=" + _stddev->str()
           + ")";
}
}  // namespace lsst::gauss2d::fit
