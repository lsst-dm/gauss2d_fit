//
// Created by dtaranu on 4/8/21.
//

#ifndef GAUSS2D_FIT_PARAMETERS_H
#define GAUSS2D_FIT_PARAMETERS_H

#include <cfloat>
#include <iostream>

#include "parameters/parameter.h"
#include "parameters/unit.h"

namespace gauss2d
{
namespace fit
{

class UnitNone : public parameters::Unit {
public:
    std::string get_name() const {
        static const std::string name_none = "None";
        return name_none; 
    }
};

static const UnitNone unit_none {};

template<typename T>
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
