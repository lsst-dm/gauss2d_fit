#ifndef GAUSS2DFIT_PARAMETRICGAUSSIAN1D_H
#define GAUSS2DFIT_PARAMETRICGAUSSIAN1D_H

#include <map>

#include "gauss2d/object.h"

#include "math.h"
#include "param_defs.h"
#include "parameters.h"

namespace gauss2d::fit {

class ParametricGaussian1D : public Object {
private:
    std::shared_ptr<MeanParameter> _mean;
    std::shared_ptr<StdDevParameter> _stddev;

public:
    double get_mean() const;
    double get_stddev() const;

    MeanParameter& get_mean_parameter() const;
    StdDevParameter& get_stddev_parameter() const;

    void set_mean(double value);
    void set_stddev(double value);

    std::string repr(bool name_keywords = false) const override;
    std::string str() const override;

    ParametricGaussian1D(std::shared_ptr<MeanParameter> mean = nullptr,
                         std::shared_ptr<StdDevParameter> stddev = nullptr);

    ~ParametricGaussian1D();
};

}  // namespace gauss2d::fit

#endif  // GAUSS2DFIT_PARAMETRICGAUSSIAN1D_H