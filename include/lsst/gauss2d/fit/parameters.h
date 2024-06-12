//
// Created by dtaranu on 4/8/21.
//

#ifndef LSST_GAUSS2D_FIT_PARAMETERS_H
#define LSST_GAUSS2D_FIT_PARAMETERS_H

#include <cfloat>
#include <iostream>

#include "lsst/modelfit/parameters.h"

namespace parameters = lsst::modelfit::parameters;

namespace lsst::gauss2d::fit {

class UnitNone : public parameters::Unit {
public:
    std::string get_name() const override {
        static const std::string name_none = "None";
        return name_none;
    }
};

static const UnitNone unit_none{};

/// A Parameter representing a size (i.e. a physical length)
struct SizeParameterD {
    virtual double get_size() const = 0;
    virtual void set_size(double size) = 0;

    virtual ~SizeParameterD() = default;
};

struct SizeXParameterD : public SizeParameterD {
    virtual ~SizeXParameterD() = default;
};

struct SizeYParameterD : public SizeParameterD {
    virtual ~SizeYParameterD() = default;
};

template <typename T>
using Param = parameters::Parameter<double, T>;

struct CentroidXParameterD : public Param<CentroidXParameterD> {
    static inline const std::string _desc = "Centroid (x)";
    static inline const std::string _name = "cen_x";
    using Param<CentroidXParameterD>::Parameter;
};

struct CentroidYParameterD : public Param<CentroidYParameterD> {
    static inline const std::string _desc = "Centroid (y)";
    static inline const std::string _name = "cen_y";
    using Param<CentroidYParameterD>::Parameter;
};

struct IntegralParameterD : public Param<IntegralParameterD> {
public:
    static inline const bool _linear = true;
    static inline constexpr double _min = 0.;
    static inline const std::string _desc = "Gaussian integral (total integrated weight/flux)";
    static inline const std::string _name = "integral";
    using Param<IntegralParameterD>::Parameter;
};

struct MeanParameterD : public Param<MeanParameterD> {
    static inline const std::string _desc = "Gaussian (1D normal) mean";
    static inline const std::string _name = "mean";
    using Param<MeanParameterD>::Parameter;
};

struct MoffatConcentrationParameterD : public Param<MoffatConcentrationParameterD> {
    static inline constexpr double _min = 1.;
    static inline constexpr double _default = 2.5;
    static inline const std::string _desc = "Moffat concentration (beta)";
    static inline const std::string _name = "con_moffat";
    using Param<MoffatConcentrationParameterD>::Parameter;
};

struct ProperFractionParameterD : public Param<ProperFractionParameterD> {
public:
    static inline constexpr double _min = 0.;
    static inline constexpr double _max = 1.;
    static inline const std::string _desc = "Proper fraction (0 <= x <= 1)";
    static inline const std::string _name = "proper_fraction";
    using Param<ProperFractionParameterD>::Parameter;
};

/// A generic scale radius, for profiles without specific names like "effective radius"
struct RadiusScaleParameterD : public Param<RadiusScaleParameterD> {
    static inline constexpr double _min = 0.;
    static inline constexpr double _default = 1.;
    static inline const std::string _desc = "Scale radius";
    static inline const std::string _name = "r_scale";
    using Param<RadiusScaleParameterD>::Parameter;
};

struct ReffXParameterD : public Param<ReffXParameterD>, SizeXParameterD {
    double get_size() const override { return this->get_value(); }
    void set_size(double size) override { this->set_value(size); }

    static inline constexpr double _min = 0.;
    static inline const std::string _desc = "Sersic effective radius (x)";
    static inline const std::string _name = "reff_x";
    using Param<ReffXParameterD>::Parameter;
};

struct ReffYParameterD : public Param<ReffYParameterD>, SizeYParameterD {
    double get_size() const override { return this->get_value(); }
    void set_size(double size) override { this->set_value(size); }

    static inline constexpr double _min = 0.;
    static inline const std::string _desc = "Sersic effective radius (y)";
    static inline const std::string _name = "reff_y";
    using Param<ReffYParameterD>::Parameter;
};

struct RhoParameterD : public Param<RhoParameterD> {
    static inline constexpr double _min = -1.;
    static inline constexpr double _default = 0.;
    static inline constexpr double _max = 1.;
    static inline const std::string _desc = "Gaussian correlation (rho)";
    static inline const std::string _name = "rho";
    using Param<RhoParameterD>::Parameter;
};

struct SersicIndexParameterD : public Param<SersicIndexParameterD> {
    static inline constexpr double _min = 0.;
    static inline constexpr double _default = 0.5;
    static inline const std::string _desc = "Sersic index";
    static inline const std::string _name = "n_ser";
    using Param<SersicIndexParameterD>::Parameter;
};

struct SigmaXParameterD : public Param<SigmaXParameterD>, SizeXParameterD {
    double get_size() const override { return this->get_value(); }
    void set_size(double size) override { this->set_value(size); }

    static inline constexpr double _min = 0.;
    static inline const std::string _desc = "Gaussian sigma (x)";
    static inline const std::string _name = "sigma_x";
    using Param<SigmaXParameterD>::Parameter;
};

struct SigmaYParameterD : public Param<SigmaYParameterD>, SizeYParameterD {
    double get_size() const override { return this->get_value(); }
    void set_size(double size) override { this->set_value(size); }

    static inline constexpr double _min = 0.;
    static inline const std::string _desc = "Gaussian sigma (y)";
    static inline const std::string _name = "sigma_y";
    using Param<SigmaYParameterD>::Parameter;
};

struct StdDevParameterD : public Param<StdDevParameterD> {
    static inline constexpr double _min = 0.;
    static inline constexpr double _default = 1.0;
    static inline const std::string _desc = "Gaussian (1D normal) standard deviation (sigma)";
    static inline const std::string _name = "stddev";
    using Param<StdDevParameterD>::Parameter;
};

}  // namespace lsst::gauss2d::fit

#endif  // LSST_GAUSS2D_FIT_PARAMETERS_H
