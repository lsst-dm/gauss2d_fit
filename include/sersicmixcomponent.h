#ifndef GAUSS2D_FIT_SERSICMIXCOMPONENT_H
#define GAUSS2D_FIT_SERSICMIXCOMPONENT_H

#include <parameters/parameter.h>

#include "channel.h"
#include "ellipticalcomponent.h"
#include "integralmodel.h"
#include "param_defs.h"
#include "param_filter.h"
#include "sersicmix.h"
#include "sersicparametricellipse.h"

namespace gauss2d::fit
{

class SersicMixComponentIndexParameter : public SersicIndexParameter {
private:
    std::vector<IntegralSize> _integralsizes;
    std::vector<IntegralSize> _integralsizes_derivs;
    const std::shared_ptr<const SersicMixInterpolator> _interpolator;

    void _set_ratios(double sersicindex);

public:
    double get_integralratio(unsigned short index) const;
    double get_integralratio_deriv(unsigned short index) const;
    const parameters::Limits<double> & get_limits_maximal() const override;
    double get_min() const override { return 0.5; }
    double get_max() const override { return 8.0; }
    double get_sizeratio(unsigned short index) const;
    double get_sizeratio_deriv(unsigned short index) const;

    unsigned short order;

    void set_value(double value) override;
    void set_value_transformed(double value_transformed) override;

    explicit SersicMixComponentIndexParameter(
        double value = _get_default(),
        std::shared_ptr<const parameters::Limits<double>> limits = nullptr,
        std::shared_ptr<const parameters::Transform<double>> transform = nullptr,
        std::shared_ptr<const parameters::Unit> unit = nullptr,
        bool fixed = false,
        std::string label = "",
        std::shared_ptr<const SersicMixInterpolator> interpolator = nullptr
    );
};

// TODO: Revisit the necessity of this class
class SersicParametricEllipseHolder {
public:
    std::shared_ptr<SersicParametricEllipse> _ellipsedata;

    explicit SersicParametricEllipseHolder(std::shared_ptr<SersicParametricEllipse> ellipse = nullptr)
    : _ellipsedata(std::move(ellipse)) {
        if(_ellipsedata == nullptr) _ellipsedata = std::make_shared<SersicParametricEllipse>();
    }
};

class SersicMixComponent : private SersicParametricEllipseHolder, public EllipticalComponent {
private:
    class SersicMixGaussianComponent;
    std::shared_ptr<SersicMixComponentIndexParameter> _sersicindex;
    std::map<std::reference_wrapper<const Channel>,
        std::vector<std::unique_ptr<SersicMixGaussianComponent>>> _gaussians;

public:
    void add_extra_param_map(const Channel & channel, ExtraParamMap & map_extra, const GradParamMap & map_grad, ParameterMap & offsets
        ) const override;
    void add_extra_param_factors(const Channel & channel, ExtraParamFactors & factors) const override;
    void add_grad_param_map(const Channel & channel, GradParamMap & map, ParameterMap & offsets
        ) const override;
    void add_grad_param_factors(const Channel & channel, GradParamFactors & factor) const override;
    
    std::unique_ptr<const gauss2d::Gaussians> get_gaussians(const Channel & channel) const override;
    size_t get_n_gaussians(const Channel & channel) const override;

    ParamRefs & get_parameters(ParamRefs & params, ParamFilter * filter = nullptr) const override;
    ParamCRefs & get_parameters_const(ParamCRefs & params, ParamFilter * filter = nullptr) const override;

    static const size_t N_PARAMS = N_PARAMS_GAUSS2D + 1;

    void set_extra_param_factors(const Channel & channel, ExtraParamFactors & factors, size_t index) const override;
    void set_grad_param_factors(const Channel & channel, GradParamFactors & factors, size_t index) const override;

    std::string repr(bool name_keywords = false) const override;
    std::string str() const override;
    
    explicit SersicMixComponent(
        std::shared_ptr<SersicParametricEllipse> ellipse = nullptr,
        std::shared_ptr<CentroidParameters> centroid = nullptr,
        std::shared_ptr<IntegralModel> integralmodel = nullptr,
        std::shared_ptr<SersicMixComponentIndexParameter> sersicindex = nullptr
    );
    ~SersicMixComponent();
};
} // namespace gauss2d::fit

#endif
