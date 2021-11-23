//
// Created by dtaranu on 4/8/21.
//

#ifndef GAUSS2D_FIT_PARAMETERS_H
#define GAUSS2D_FIT_PARAMETERS_H

#include <cfloat>
#include <iostream>

#ifndef PARAMETERS_PARAMETER_H
#include "parameters/parameter.h"
#endif

namespace gauss2d
{
namespace fit
{
class Unit {
public:
    virtual std::string name() const = 0;
    virtual ~Unit() = default;
};

static const std::string name_unit_none = "None";
class UnitNone : public Unit {
public:
    std::string name() const { return name_unit_none; }
};

static const UnitNone unit_none {};

template<typename T>
using Param = parameters::Parameter<double, T, Unit>;

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
private:
    bool _is_ratio = false;

public:
    bool get_is_ratio() const { return _is_ratio; }
    void set_is_ratio(bool is_ratio) { _is_ratio = is_ratio; }

    static inline constexpr double _min = 0.;
    static inline const std::string _desc = "Gaussian integral (total integrated weight/flux)";
    static inline const std::string _name = "integral";
    using Param<IntegralParameter>::Parameter;
};

struct MoffatConcentrationParameter : public Param<MoffatConcentrationParameter> {
    static inline constexpr double _min = 1.;
    static inline constexpr double _default = 2.5;
    static inline const std::string _desc = "Moffat concentration (beta)";
    static inline const std::string _name = "con_moffat";
    using Param<MoffatConcentrationParameter>::Parameter;
};

struct RadiusScaleParameter : public Param<RadiusScaleParameter> {
    static inline constexpr double _min = 0.;
    static inline constexpr double _default = 1.;
    static inline const std::string _desc = "Scale radius";
    static inline const std::string _name = "r_scale";
    using Param<RadiusScaleParameter>::Parameter;
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

struct SigmaXParameter : public Param<SigmaXParameter> {
    static inline constexpr double _min = 0.;
    static inline const std::string _desc = "Gaussian sigma (x)";
    static inline const std::string _name = "sigma_x";
    using Param<SigmaXParameter>::Parameter;
};

struct SigmaYParameter : public Param<SigmaYParameter> {
    static inline constexpr double _min = 0.;
    static inline const std::string _desc = "Gaussian sigma (y)";
    static inline const std::string _name = "sigma_y";
    using Param<SigmaYParameter>::Parameter;
};

}
}

#endif //GAUSS2DFIT_PARAMETERS_H
