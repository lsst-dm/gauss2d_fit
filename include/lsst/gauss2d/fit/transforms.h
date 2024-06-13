#ifndef LSST_GAUSS2D_FIT_TRANSFORMS_H
#define LSST_GAUSS2D_FIT_TRANSFORMS_H

#include <cmath>
#include <iostream>
#include <memory>

#include "lsst/modelfit/parameters.h"

namespace lsst::gauss2d::fit {
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wnon-virtual-dtor"

namespace parameters = lsst::modelfit::parameters;
typedef parameters::Transform<double> Transform;

struct InverseTransform : public Transform {
    std::string description() const override { return "Inverse transform"; }
    std::string repr(bool name_keywords = false,
                     const std::string_view& namespace_separator
                     = parameters::Object::CC_NAMESPACE_SEPARATOR) const override {
        return parameters::type_name_str<InverseTransform>(false, namespace_separator) + "()";
    }
    std::string str() const override { return parameters::type_name_str<InverseTransform>(true) + "()"; }

    inline double derivative(double x) const override { return 1 / (x * x); }
    inline double forward(double x) const override { return 1 / x; }
    inline double reverse(double x) const override { return 1 / x; }
};

struct JanskyToABMagTransform : public Transform {
    static inline const double f_nu_0 = 3630.780547701002879554236770479;

    std::string description() const override { return "jansky to AB magnitude transform"; }
    std::string repr(bool name_keywords = false,
                     const std::string_view& namespace_separator
                     = parameters::Object::CC_NAMESPACE_SEPARATOR) const override {
        return parameters::type_name_str<JanskyToABMagTransform>(false, namespace_separator) + "()";
    }
    std::string str() const override {
        return parameters::type_name_str<JanskyToABMagTransform>(true) + "()";
    }

    inline double derivative(double x) const override {
        return -1.08573620475812959718098227313021197915 / x;
    }
    inline double forward(double x) const override { return -2.5 * log10(x / f_nu_0); }
    inline double reverse(double x) const override { return f_nu_0 * pow10(-0.4 * x); }
};

struct NanojanskyToABMagTransform : public JanskyToABMagTransform {
    std::string description() const override { return "nanojansky to AB magnitude transform"; }
    std::string repr(bool name_keywords = false,
                     const std::string_view& namespace_separator
                     = parameters::Object::CC_NAMESPACE_SEPARATOR) const override {
        return parameters::type_name_str<NanojanskyToABMagTransform>(false, namespace_separator) + "()";
    }
    std::string str() const override {
        return parameters::type_name_str<NanojanskyToABMagTransform>(true) + "()";
    }

    inline double derivative(double x) const override { return JanskyToABMagTransform::derivative(x); }
    inline double forward(double x) const override { return JanskyToABMagTransform::forward(x * 1e-9); }
    inline double reverse(double x) const override { return 1e9 * JanskyToABMagTransform::reverse(x); }
};

struct LogTransform : public Transform {
    std::string description() const override { return "Natural (base e) logarithmic transform"; }
    std::string repr(bool name_keywords = false,
                     const std::string_view& namespace_separator
                     = parameters::Object::CC_NAMESPACE_SEPARATOR) const override {
        return parameters::type_name_str<LogTransform>(false, namespace_separator) + "()";
    }
    std::string str() const override { return parameters::type_name_str<LogTransform>(true) + "()"; }

    inline double derivative(double x) const override { return 1 / x; }
    inline double forward(double x) const override { return log(x); }
    inline double reverse(double x) const override { return exp(x); }
};

struct Log10Transform : public Transform {
    std::string description() const override { return "Base 10 logarithmic transform"; }
    std::string repr(bool name_keywords = false,
                     const std::string_view& namespace_separator
                     = parameters::Object::CC_NAMESPACE_SEPARATOR) const override {
        return parameters::type_name_str<Log10Transform>(false, namespace_separator) + "()";
    }
    std::string str() const override { return parameters::type_name_str<Log10Transform>(true) + "()"; }

    inline double derivative(double x) const override {
        return 0.434294481903251827651128918916605082294397 / x;
    }
    inline double forward(double x) const override { return log10(x); }
    inline double reverse(double x) const override { return pow10(x); }
};

struct LogitTransform : public Transform {
    std::string description() const override { return "Logit transform"; }
    std::string repr(bool name_keywords = false,
                     const std::string_view& namespace_separator
                     = parameters::Object::CC_NAMESPACE_SEPARATOR) const override {
        return parameters::type_name_str<LogitTransform>(false, namespace_separator) + "()";
    }
    std::string str() const override { return parameters::type_name_str<LogitTransform>(true) + "()"; }

    inline double derivative(double x) const override { return 1 / x + 1 / (1 - x); }
    inline double forward(double x) const override { return log(x / (1 - x)); }
    inline double reverse(double x) const override { return 1 / (1 + exp(-x)); }
};

class LogitLimitedTransform : public Transform {
public:
    explicit LogitLimitedTransform(std::shared_ptr<parameters::Limits<double>> limits, double factor = 1) {
        set_limits(std::move(limits));
        set_factor(factor);
        _set_range();
    }

    std::string description() const override { return "Logit limited (to finite range) transform"; }
    std::string repr(bool name_keywords = false,
                     const std::string_view& namespace_separator
                     = parameters::Object::CC_NAMESPACE_SEPARATOR) const override {
        return parameters::type_name_str<LogitLimitedTransform>(false, namespace_separator) + "("
               + (name_keywords ? "limits=" : "") + _limits->repr(name_keywords, namespace_separator) + ", "
               + (name_keywords ? "factor=" : "") + std::to_string(_factor) + ")";
    }
    std::string str() const override {
        return parameters::type_name_str<LogitLimitedTransform>(true) + "(limits=" + _limits->str()
               + ", factor=" + std::to_string(_factor) + ")";
    }

    double get_factor() const { return _factor; }
    parameters::Limits<double>& get_limits() const { return *_limits; }

    double derivative(double x) const override;
    double forward(double x) const override;
    double reverse(double x) const override;

    void set_factor(double factor) {
        if (!(factor > 0))
            throw std::invalid_argument("LogitLimitedTransform factor=" + std::to_string(factor) + " !>0");
        _factor = factor;
    }

    void set_limits(std::shared_ptr<parameters::Limits<double>> limits) {
        _limits = (limits == nullptr) ? std::make_shared<parameters::Limits<double>>() : std::move(limits);
        _set_range();
    }

private:
    std::shared_ptr<parameters::Limits<double>> _limits;
    double _factor;
    double _range;

    inline void _set_range() { _range = _limits->get_max() - _limits->get_min(); }
};

template <class T>
std::shared_ptr<T> get_transform_default() {
    static T transform_default{};
    static std::shared_ptr<T> ptr{std::shared_ptr<T>{}, &transform_default};
    return ptr;
}
#pragma GCC diagnostic pop
}  // namespace lsst::gauss2d::fit

#endif  // LSST_GAUSS2D_FIT_TRANSFORMS_H
