//
// Created by dtaranu on 4/8/21.
//

#ifndef GAUSS2D_FIT_TRANSFORMS_H
#define GAUSS2D_FIT_TRANSFORMS_H

#include <cmath>
#include <iostream>
#include <memory>

#ifndef PARAMETERS_LIMITS_H
#include "parameters/limits.h"
#endif

#ifndef PARAMETERS_TRANSFORM_H
#include "parameters/transform.h"
#endif

namespace gauss2d
{
namespace fit
{

typedef parameters::Transform<double> Transform;

struct InverseTransform : public Transform
{
    std::string description() const { return "Inverse transform"; }
    std::string str() const { return "InverseTransform"; }

    inline double derivative(double x) const { return 1/(x*x); }
    inline double forward(double x) const { return 1/x; }
    inline double reverse(double x) const { return 1/x; }
};

struct LogTransform : public Transform
{
    std::string description() const { return "Natural (base e) logarithmic transform"; }
    std::string str() const { return "LogTransform"; }

    inline double derivative(double x) const { return 1/x; }
    inline double forward(double x) const { return log(x); }
    inline double reverse(double x) const { return exp(x); }
};

struct Log10Transform : public Transform
{
    std::string description() const { return "Base 10 logarithmic transform"; }
    std::string str() const { return "Log10Transform"; }

    inline double derivative(double x) const { return 0.434294481903251827651128918916605082294397/x; }
    inline double forward(double x) const { return log10(x); }
    inline double reverse(double x) const { return pow10(x); }
};

struct LogitTransform : public Transform
{
    std::string description() const { return "Logit transform"; }
    std::string str() const { return "LogitTransform"; }

    inline double derivative(double x) const { return 1/x + 1/(1 - x); }
    inline double forward(double x) const { return log(x/(1 - x)); }
    inline double reverse(double x) const { return 1/(1 + exp(-x)); }
};

class LogitLimitedTransform : public Transform
{
private:
    std::shared_ptr<parameters::Limits<double>> _limits;
    double _factor;
    double _range;

    inline void _set_range() { _range = _limits->get_max() - _limits->get_min(); }

public:
    std::string description() const { return "Logit limited (to finite range) transform"; }
    std::string str() const {
        return std::string(parameters::type_name<LogitLimitedTransform>()) + "(limits=" + _limits->str() 
            + ", factor=" + std::to_string(_factor) + ")";
    }

    double get_factor() const { return _factor; }
    parameters::Limits<double> & get_limits() const { return *_limits; }

    double derivative(double x) const;
    double forward(double x) const;
    double reverse(double x) const;

    void set_factor(double factor) {
        if(!(factor > 0)) throw std::invalid_argument("LogitLimitedTransform factor="
            + std::to_string(factor) + " !>0");
        _factor = factor;
    }

    void set_limits(std::shared_ptr<parameters::Limits<double>> limits) {
        _limits = (limits == nullptr) ? std::make_shared<parameters::Limits<double>>() : std::move(limits);
        _set_range();
    }

    LogitLimitedTransform(
        std::shared_ptr<parameters::Limits<double>> limits,
        double factor = 1
    ) {
        set_limits(limits);
        set_factor(factor);
        _set_range();
    }
};
}
}

#endif //GAUSS2DFIT_TRANSFORMS_H
