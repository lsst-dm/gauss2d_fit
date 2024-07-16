#ifndef LSST_GAUSS2D_FIT_PARAMETRICGAUSSIAN1D_H
#define LSST_GAUSS2D_FIT_PARAMETRICGAUSSIAN1D_H

#include <map>

#include "lsst/gauss2d/object.h"

#include "math.h"
#include "param_defs.h"
#include "parameters.h"

namespace lsst::gauss2d::fit {

class ParametricGaussian1D : public Object {
public:
    explicit ParametricGaussian1D(std::shared_ptr<MeanParameterD> mean = nullptr,
                                  std::shared_ptr<StdDevParameterD> stddev = nullptr);

    ~ParametricGaussian1D();

    double get_mean() const;
    double get_stddev() const;

    MeanParameterD& get_mean_parameter() const;
    StdDevParameterD& get_stddev_parameter() const;

    void set_mean(double value);
    void set_stddev(double value);

    std::string repr(bool name_keywords = false,
                     std::string_view namespace_separator = Object::CC_NAMESPACE_SEPARATOR) const override;
    std::string str() const override;

private:
    std::shared_ptr<MeanParameterD> _mean;
    std::shared_ptr<StdDevParameterD> _stddev;
};

}  // namespace lsst::gauss2d::fit

#endif  // GAUSS2D_FIT_PARAMETRICGAUSSIAN1D_H