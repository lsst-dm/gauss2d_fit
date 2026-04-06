#ifndef LSST_GAUSS2D_FIT_TRANSFORMS_H
#define LSST_GAUSS2D_FIT_TRANSFORMS_H

#include <cmath>
#include <iostream>
#include <memory>

#include "lsst/gauss2d/string_utils.h"

#include "lsst/modelfit/parameters.h"

#include "util.h"

namespace lsst::gauss2d::fit {
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wnon-virtual-dtor"

namespace parameters = lsst::modelfit::parameters;

template <typename C, typename T>
std::string stripped_template_type_name_str(
    std::string prefix = "",
    bool strip_namespace = false,
    const std::string_view& namespace_str = "::") {
    std::string name = parameters::type_name_str<C>(strip_namespace, namespace_str);
    std::string to_replace = prefix + "<" + parameters::type_name_str<T>() + ">";
    name = replace_all_none(name, to_replace);
    return name;
}

template <typename T>
struct InverseTransform_ : public parameters::Transform<T> {
    std::string description() const override { return "Inverse transform"; }
    std::string repr(bool name_keywords = false,
                     const std::string_view& namespace_separator
                     = parameters::Object::CC_NAMESPACE_SEPARATOR) const override {
        return stripped_template_type_name_str<InverseTransform_<T>, T>("_", false, namespace_separator) + suffix_type_str<T>() + "()";
    }
    std::string str() const override { return stripped_template_type_name_str<InverseTransform_<T>, T>("_", true) + suffix_type_str<T>() + "()"; }

    inline T derivative(T x) const override { return 1 / (x * x); }
    inline T forward(T x) const override { return 1 / x; }
    inline T reverse(T x) const override { return 1 / x; }
};

template <typename T>
struct JanskyToABMagTransform_ : public parameters::Transform<T> {
    static inline const double f_nu_0 = 3630.780547701002879554236770479;

    std::string description() const override { return "jansky to AB magnitude transform"; }
    std::string repr(bool name_keywords = false,
                     const std::string_view& namespace_separator
                     = parameters::Object::CC_NAMESPACE_SEPARATOR) const override {
        return stripped_template_type_name_str<JanskyToABMagTransform_<T>, T>("_", false, namespace_separator) + suffix_type_str<T>() + "()";
    }
    std::string str() const override {
        return stripped_template_type_name_str<JanskyToABMagTransform_<T>, T>("_", true) + suffix_type_str<T>() + "()";
    }

    inline T derivative(T x) const override {
        return -1.08573620475812959718098227313021197915 / x;
    }
    inline T forward(T x) const override { return -2.5 * log10(x / f_nu_0); }
    inline T reverse(T x) const override { return f_nu_0 * pow(10.0, -0.4 * x); }
};

template <typename T>
struct NanojanskyToABMagTransform_ : public JanskyToABMagTransform_<T> {
    std::string description() const override { return "nanojansky to AB magnitude transform"; }
    std::string repr(bool name_keywords = false,
                     const std::string_view& namespace_separator
                     = parameters::Object::CC_NAMESPACE_SEPARATOR) const override {
        return stripped_template_type_name_str<NanojanskyToABMagTransform_<T>, T>("_", false, namespace_separator) + suffix_type_str<T>() + "()";
    }
    std::string str() const override {
        return stripped_template_type_name_str<NanojanskyToABMagTransform_<T>, T>("_", true) + suffix_type_str<T>() + "()";
    }

    inline T derivative(T x) const override { return JanskyToABMagTransform_<T>::derivative(x); }
    inline T forward(T x) const override { return JanskyToABMagTransform_<T>::forward(x * 1e-9); }
    inline T reverse(T x) const override { return 1e9 * JanskyToABMagTransform_<T>::reverse(x); }
};

template <typename T>
struct LogTransform_ : public parameters::Transform<T> {
    std::string description() const override { return "Natural (base e) logarithmic transform"; }
    std::string repr(bool name_keywords = false,
                     const std::string_view& namespace_separator
                     = parameters::Object::CC_NAMESPACE_SEPARATOR) const override {
        return stripped_template_type_name_str<LogTransform_<T>, T>("_", false, namespace_separator) + suffix_type_str<T>() + "()";
    }
    std::string str() const override { return stripped_template_type_name_str<LogTransform_<T>, T>("_", true) + suffix_type_str<T>() + "()"; }

    inline T derivative(T x) const override { return 1 / x; }
    inline T forward(T x) const override { return log(x); }
    inline T reverse(T x) const override { return exp(x); }
};

template <typename T>
struct Log10Transform_ : public parameters::Transform<T> {
    std::string description() const override { return "Base 10 logarithmic transform"; }
    std::string repr(bool name_keywords = false,
                     const std::string_view& namespace_separator
                     = parameters::Object::CC_NAMESPACE_SEPARATOR) const override {
        return stripped_template_type_name_str<Log10Transform_<T>, T>("_", false, namespace_separator) + suffix_type_str<T>() + "()";
    }
    std::string str() const override { return stripped_template_type_name_str<Log10Transform_<T>, T>("_", true) + suffix_type_str<T>() + "()"; }

    inline T derivative(T x) const override {
        return 0.434294481903251827651128918916605082294397 / x;
    }
    inline T forward(T x) const override { return log10(x); }
    inline T reverse(T x) const override { return pow(10., x); }
};

template <typename T>
struct LogitTransform_ : public parameters::Transform<T> {
    std::string description() const override { return "Logit transform"; }
    std::string repr(bool name_keywords = false,
                     const std::string_view& namespace_separator
                     = parameters::Object::CC_NAMESPACE_SEPARATOR) const override {
        return stripped_template_type_name_str<LogitTransform_<T>, T>("_", false, namespace_separator) + suffix_type_str<T>() + "()";
    }
    std::string str() const override { return stripped_template_type_name_str<LogitTransform_<T>, T>("_", true) + suffix_type_str<T>() + "()"; }

    inline T derivative(T x) const override { return 1 / x + 1 / (1 - x); }
    inline T forward(T x) const override { return log(x / (1 - x)); }
    inline T reverse(T x) const override { return 1 / (1 + exp(-x)); }
};

template <typename T>
class LogitLimitedTransform_ : public parameters::Transform<T> {
public:
    explicit LogitLimitedTransform_(std::shared_ptr<parameters::Limits<T>> limits, double factor = 1) {
        set_limits(std::move(limits));
        set_factor(factor);
        _set_range();
    }

    std::string description() const override { return "Logit limited (to finite range) transform"; }
    std::string repr(bool name_keywords = false,
                     const std::string_view& namespace_separator
                     = parameters::Object::CC_NAMESPACE_SEPARATOR) const override {
        return stripped_template_type_name_str<LogitLimitedTransform_<T>, T>("_", false, namespace_separator) + "("
               + (name_keywords ? "limits=" : "") + _limits->repr(name_keywords, namespace_separator) + ", "
               + (name_keywords ? "factor=" : "") + std::to_string(_factor) + ")";
    }
    std::string str() const override {
        return stripped_template_type_name_str<LogitLimitedTransform_<T>, T>("_", true) + "(limits=" + _limits->str()
               + ", factor=" + std::to_string(_factor) + ")";
    }

    T get_factor() const { return _factor; }
    parameters::Limits<T>& get_limits() const { return *_limits; }

    T derivative(T x) const override {
        T y = (x - _limits->get_min()) / _range;
        if (y == 1) {
            return std::numeric_limits<T>::infinity();;
        } else if (y == 0) {
            return -std::numeric_limits<T>::infinity();;
        }
        return (1 / y + 1 / (1 - y)) * _factor / _range;
    }

    T forward(T x) const override {
        T min = _limits->get_min();
        T max = _limits->get_max();
        if (x == min)
            return -std::numeric_limits<T>::infinity();
        else if (x == max)
            return std::numeric_limits<T>::infinity();
        double y = (x - min) / _range;
        if (!(y < 1) || !(y > 0)) return nan("");
        return log(y / (1 - y)) * _factor;
    }

    T reverse(T x) const override {
        T y = -x * _factor;
        // both exp(y) and 1/y could blow up near +/- inf
        if (y > _max_good()) return _limits->get_min();
        y = 1 + ((y < -_max_good()) ? 0 : exp(y));
        return _range / y + _limits->get_min();
    }

    void set_factor(double factor) {
        if (!(factor > 0))
            throw std::invalid_argument("LogitLimitedTransform" + suffix_type_str<T>() + " factor=" + std::to_string(factor) + " !>0");
        _factor = factor;
    }

    void set_limits(std::shared_ptr<parameters::Limits<T>> limits) {
        _limits = (limits == nullptr) ? std::make_shared<parameters::Limits<T>>() : std::move(limits);
        _set_range();
    }

private:
    std::shared_ptr<parameters::Limits<T>> _limits;
    T _factor;
    T _range;

    static const T _max_good() {
        return log(0.999) + log(std::numeric_limits<T>::max());
    }

    inline void _set_range() { _range = _limits->get_max() - _limits->get_min(); }
};

typedef InverseTransform_<double> InverseTransformD;
[[deprecated("Use InverseTransformD instead")]]
typedef InverseTransform_<double> InverseTransform;

typedef JanskyToABMagTransform_<double> JanskyToABMagTransformD;
[[deprecated("Use JanskyToABMagTransformD instead")]]
typedef JanskyToABMagTransform_<double> JanskyToABMagTransform;

typedef NanojanskyToABMagTransform_<double> NanojanskyToABMagTransformD;
[[deprecated("Use NanojanskyToABMagTransformD instead")]]
typedef NanojanskyToABMagTransform_<double> NanojanskyToABMagTransform;

typedef LogTransform_<double> LogTransformD;
[[deprecated("Use LogTransformD instead")]]
typedef LogTransform_<double> LogTransform;

typedef Log10Transform_<double> Log10TransformD;
[[deprecated("Use Log10TransformD instead")]]
typedef Log10Transform_<double> Log10Transform;

typedef LogitTransform_<double> LogitTransformD;
[[deprecated("Use LogitTransformD instead")]]
typedef LogitTransform_<double> LogitTransform;

typedef LogitLimitedTransform_<double> LogitLimitedTransformD;
[[deprecated("Use LogitLimitedTransformD instead")]]
typedef LogitLimitedTransform_<double> LogitLimitedTransform;

template <class T>
std::shared_ptr<T> get_transform_default() {
    static T transform_default{};
    static std::shared_ptr<T> ptr{std::shared_ptr<T>{}, &transform_default};
    return ptr;
}
#pragma GCC diagnostic pop
}  // namespace lsst::gauss2d::fit

#endif  // LSST_GAUSS2D_FIT_TRANSFORMS_H
