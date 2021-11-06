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
template<typename T>
using Param = parameters::Parameter<double, T>;

struct IntegralParameter : public Param<IntegralParameter> {
    static inline constexpr double _min = 0.;
    static inline const std::string _desc = "Gaussian integral (total integrated weight/flux)";
    static inline const std::string _name = "integral";
    using Param<IntegralParameter>::Parameter;
};

struct RhoParameter : public Param<RhoParameter> {
    static inline constexpr double _min = 0.;
    static inline constexpr double _default = 1.;
    static inline constexpr double _max = 1.;
    static inline const std::string _desc = "Gaussian correlation (rho)";
    static inline const std::string _name = "rho";
    using Param<RhoParameter>::Parameter;
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

}
}

#endif //GAUSS2DFIT_PARAMETERS_H
