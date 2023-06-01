#include "parameters.h"
#include "parametricgaussian1d.h"

namespace gauss2d::fit {

double ParametricGaussian1D::get_mean() const { return this->_mean->get_value(); }

double ParametricGaussian1D::get_stddev() const { return this->_stddev->get_value(); }

const MeanParameter &ParametricGaussian1D::get_mean_parameter() const { return *this->_mean; }

const StdDevParameter &ParametricGaussian1D::get_stddev_parameter() const { return *this->_stddev; }

void ParametricGaussian1D::set_mean(double value) { this->_mean->set_value(value); }

void ParametricGaussian1D::set_stddev(double value) { this->_stddev->set_value(value); }

std::string ParametricGaussian1D::repr(bool name_keywords) const {
    return std::string("ParametricGaussian1D(") + (name_keywords ? "mean=" : "") + _mean->repr(name_keywords)
           + (name_keywords ? ", stddev=" : "") + _stddev->repr(name_keywords) + ")";
}

std::string ParametricGaussian1D::str() const {
    return std::string("ParametricGaussian1D(mean=") + _mean->str() + ", stddev=" + _stddev->str() + ")";
}

ParametricGaussian1D::ParametricGaussian1D(std::shared_ptr<MeanParameter> mean,
                                           std::shared_ptr<StdDevParameter> stddev)
        : _mean(std::move(mean)), _stddev(std::move(stddev)) {
    if (_mean == nullptr) _mean = std::make_shared<MeanParameter>();
    if (_stddev == nullptr) _stddev = std::make_shared<StdDevParameter>();
}

ParametricGaussian1D::~ParametricGaussian1D() {}

}  // namespace gauss2d::fit
