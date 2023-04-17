#ifndef GAUSS2D_FIT_ELLIPTICALCOMPONENT_H
#define GAUSS2D_FIT_ELLIPTICALCOMPONENT_H

#include "centroidparameters.h"
#include "component.h"
#include "param_defs.h"
#include "param_filter.h"
#include "parametricellipse.h"

namespace gauss2d::fit {

/**
 * A Component with an elliptically-symmetric intensity profile.
 */
class EllipticalComponent : public Component {
protected:
    std::shared_ptr<ParametricEllipse> _ellipse;
    std::shared_ptr<CentroidParameters> _centroid;
    std::shared_ptr<IntegralModel> _integralmodel;

public:
    /// Get the centroid Parameter container
    const CentroidParameters& get_centroid() const;
    /// Get the ellipse Parameter container
    const ParametricEllipse& get_ellipse() const;
    const IntegralModel& get_integralmodel() const override;

    ParamRefs& get_parameters(ParamRefs& params, ParamFilter* filter = nullptr) const override;
    ParamCRefs& get_parameters_const(ParamCRefs& params, ParamFilter* filter = nullptr) const override;

    std::string repr(bool name_keywords = false) const override;
    std::string str() const override;

    /**
     * Construct an EllipticalComponent from Parameter containers
     *
     * @param ellipse The ellipse Parameter container
     * @param centroid The centroid Parameter container; default-constructed if null
     * @param integralmodel The IntegralModel; default-constructed if null
     */
    explicit EllipticalComponent(std::shared_ptr<ParametricEllipse> ellipse,
                                 std::shared_ptr<CentroidParameters> centroid = nullptr,
                                 std::shared_ptr<IntegralModel> integralmodel = nullptr);
};

}  // namespace gauss2d::fit

#endif
