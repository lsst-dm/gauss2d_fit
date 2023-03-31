#ifndef GAUSS2D_FIT_PARAMETRIC_H
#define GAUSS2D_FIT_PARAMETRIC_H

#include "gauss2d/object.h"

#include "param_defs.h"
#include "param_filter.h"

namespace gauss2d
{
namespace fit
{

class Parametric : public Object
{
public:
    virtual ParamRefs & get_parameters(ParamRefs & params, ParamFilter * filter = nullptr) const = 0;
    virtual ParamCRefs & get_parameters_const(ParamCRefs & params, ParamFilter * filter = nullptr) const = 0;
    ParamRefs get_parameters_new(ParamFilter * filter = nullptr) const {
        ParamRefs params {};
        get_parameters(params, filter);
        return params;
    }
    ParamCRefs get_parameters_const_new(ParamFilter * filter = nullptr) const {
        ParamCRefs params {};
        get_parameters_const(params, filter);
        return params;
    }

    //virtual std::string repr(bool name_keywords = false) const override = 0;
    //virtual std::string str() const override = 0;
    virtual ~Parametric() {};
};
} // namespace fit
} // namespace gauss2d

#endif
