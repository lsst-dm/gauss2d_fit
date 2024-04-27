//
// Created by dtaranu on 4/8/21.
//

#ifndef GAUSS2D_FIT_PARAMETERS_H
#define GAUSS2D_FIT_PARAMETERS_H

#include <cfloat>
#include <iostream>

#include "lsst/modelfit/parameters.h"

namespace parameters = lsst::modelfit::parameters;

namespace gauss2d::fit {

class UnitNone : public parameters::Unit {
public:
    std::string get_name() const override {
        static const std::string name_none = "None";
        return name_none;
    }
};

static const UnitNone unit_none{};

/// A Parameter representing a size (i.e. a physical length)
struct SizeParameter {
    virtual double get_size() const = 0;
    virtual void set_size(double size) = 0;

    virtual ~SizeParameter() = default;
};

struct SizeXParameter : public SizeParameter {
    virtual ~SizeXParameter() = default;
};

struct SizeYParameter : public SizeParameter {
    virtual ~SizeYParameter() = default;
};

template <typename T>
using Param = parameters::Parameter<double, T>;

struct CentroidXParameter : public Param<CentroidXParameter> {
    static inline const std::string _desc = "Centroid (x)";
    static inline const std::string _name = "cen_x";
    using Param<CentroidXParameter>::Parameter;
};

struct CentroidYParameter : public Param<CentroidYParameter> {
    static inline const std::string _desc = "Centroid (y)";
    static inline const std::string _name = "cen_y";
    using Param<CentroidYParameter>::Parameter;
};

struct IntegralParameter : public Param<IntegralParameter> {
public:
    static inline const bool _linear = true;
    static inline constexpr double _min = 0.;
    static inline const std::string _desc = "Gaussian integral (total integrated weight/flux)";
    static inline const std::string _name = "integral";
    using Param<IntegralParameter>::Parameter;
};

struct MeanParameter : public Param<MeanParameter> {
    static inline const std::string _desc = "Gaussian (1D normal) mean";
    static inline const std::string _name = "mean";
    using Param<MeanParameter>::Parameter;
};

struct MoffatConcentrationParameter : public Param<MoffatConcentrationParameter> {
    static inline constexpr double _min = 1.;
    static inline constexpr double _default = 2.5;
    static inline const std::string _desc = "Moffat concentration (beta)";
    static inline const std::string _name = "con_moffat";
    using Param<MoffatConcentrationParameter>::Parameter;
};

struct ProperFractionParameter : public Param<ProperFractionParameter> {
public:
    static inline constexpr double _min = 0.;
    static inline constexpr double _max = 1.;
    static inline const std::string _desc = "Proper fraction (0 <= x <= 1)";
    static inline const std::string _name = "proper_fraction";
    using Param<ProperFractionParameter>::Parameter;
};

/// A generic scale radius, for profiles without specific names like "effective radius"
struct RadiusScaleParameter : public Param<RadiusScaleParameter> {
    static inline constexpr double _min = 0.;
    static inline constexpr double _default = 1.;
    static inline const std::string _desc = "Scale radius";
    static inline const std::string _name = "r_scale";
    using Param<RadiusScaleParameter>::Parameter;
};

struct ReffXParameter : public Param<ReffXParameter>, SizeXParameter {
    double get_size() const override { return this->get_value(); }
    void set_size(double size) override { this->set_value(size); }

    static inline constexpr double _min = 0.;
    static inline const std::string _desc = "Sersic effective radius (x)";
    static inline const std::string _name = "reff_x";
    using Param<ReffXParameter>::Parameter;
};

struct ReffYParameter : public Param<ReffYParameter>, SizeYParameter {
    double get_size() const override { return this->get_value(); }
    void set_size(double size) override { this->set_value(size); }

    static inline constexpr double _min = 0.;
    static inline const std::string _desc = "Sersic effective radius (y)";
    static inline const std::string _name = "reff_y";
    using Param<ReffYParameter>::Parameter;
};

struct RhoParameter : public Param<RhoParameter> {
    static inline constexpr double _min = -1.;
    static inline constexpr double _default = 0.;
    static inline constexpr double _max = 1.;
    static inline const std::string _desc = "Gaussian correlation (rho)";
    static inline const std::string _name = "rho";
    using Param<RhoParameter>::Parameter;
};

struct SersicIndexParameter : public Param<SersicIndexParameter> {
    static inline constexpr double _min = 0.;
    static inline constexpr double _default = 0.5;
    static inline const std::string _desc = "Sersic index";
    static inline const std::string _name = "n_ser";
    using Param<SersicIndexParameter>::Parameter;
};

struct SigmaXParameter : public Param<SigmaXParameter>, SizeXParameter {
    double get_size() const override { return this->get_value(); }
    void set_size(double size) override { this->set_value(size); }

    static inline constexpr double _min = 0.;
    static inline const std::string _desc = "Gaussian sigma (x)";
    static inline const std::string _name = "sigma_x";
    using Param<SigmaXParameter>::Parameter;
};

struct SigmaYParameter : public Param<SigmaYParameter>, SizeYParameter {
    double get_size() const override { return this->get_value(); }
    void set_size(double size) override { this->set_value(size); }

    static inline constexpr double _min = 0.;
    static inline const std::string _desc = "Gaussian sigma (y)";
    static inline const std::string _name = "sigma_y";
    using Param<SigmaYParameter>::Parameter;
};

struct StdDevParameter : public Param<StdDevParameter> {
    static inline constexpr double _min = 0.;
    static inline constexpr double _default = 1.0;
    static inline const std::string _desc = "Gaussian (1D normal) standard deviation (sigma)";
    static inline const std::string _name = "stddev";
    using Param<StdDevParameter>::Parameter;
};

}  // namespace gauss2d::fit

#endif  // GAUSS2D_FIT_PARAMETERS_H
