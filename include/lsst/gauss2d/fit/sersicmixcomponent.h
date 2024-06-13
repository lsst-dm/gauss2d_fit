#ifndef LSST_GAUSS2D_FIT_SERSICMIXCOMPONENT_H
#define LSST_GAUSS2D_FIT_SERSICMIXCOMPONENT_H

#include "lsst/modelfit/parameters.h"

#include "channel.h"
#include "ellipticalcomponent.h"
#include "integralmodel.h"
#include "param_defs.h"
#include "param_filter.h"
#include "sersicmix.h"
#include "sersicparametricellipse.h"

namespace lsst::gauss2d::fit {
/**
 * A SersicIndexParameter for a Gaussian mixture Component.
 */
class SersicMixComponentIndexParameterD : public SersicIndexParameterD {
public:
    /// See docs for Parameter
    explicit SersicMixComponentIndexParameterD(
            double value = _get_default(), std::shared_ptr<const parameters::Limits<double>> limits = nullptr,
            std::shared_ptr<const parameters::Transform<double>> transform = nullptr,
            std::shared_ptr<const parameters::Unit> unit = nullptr, bool fixed = false,
            std::string label = "", std::shared_ptr<const SersicMixInterpolator> interpolator = nullptr);

    /// Return the integral ratio for a given Gaussian sub-component index
    double get_integralratio(unsigned short index) const;
    /// Return the integral ratio derivative for a given Gaussian sub-component index
    double get_integralratio_deriv(unsigned short index) const;
    static std::shared_ptr<const SersicMixInterpolator> get_interpolator_default(unsigned short order
                                                                                 = SERSICMIX_ORDER_DEFAULT);
    const parameters::Limits<double>& get_limits_maximal() const override;
    double get_min() const override { return 0.5; }
    double get_max() const override { return 8.0; }
    /// Return the size ratio for a given Gaussian sub-component index
    double get_sizeratio(unsigned short index) const;
    /// Return the size ratio derivative for a given Gaussian sub-component index
    double get_sizeratio_deriv(unsigned short index) const;

    InterpType get_interptype() const;
    unsigned short get_order() const;

    void set_value(double value) override;
    void set_value_transformed(double value_transformed) override;

private:
    std::vector<IntegralSize> _integralsizes;
    std::vector<IntegralSize> _integralsizes_derivs;
    const std::shared_ptr<const SersicMixInterpolator> _interpolator;

    /// Add the values for a given Sersic index
    void _set_ratios(double sersicindex);
};

// TODO: Revisit the necessity of this class
/**
 * A workaround class to store a SersicParametricEllipse. Not for direct use.
 */
class SersicParametricEllipseHolder {
public:
    /// This constructor does not need to be called by users.
    explicit SersicParametricEllipseHolder(std::shared_ptr<SersicParametricEllipse> ellipse = nullptr)
            : _ellipsedata(std::move(ellipse)) {
        if (_ellipsedata == nullptr) _ellipsedata = std::make_shared<SersicParametricEllipse>();
    }
    std::shared_ptr<SersicParametricEllipse> _ellipsedata;
};

// TODO: Add ref to derivation of weights, when published
/**
 * @brief A Gaussian mixture approximation to a Sersic profile Component.
 *
 * This Component approximates the Sersic (1968) profile, which is a widely-used
 * (in astronomy) exponentially-declining radial profile with a single shape
 * parameter. It is also functionally a generalization of the Gaussian (Sersic
 * index n=0.5), exponential (n=1) and de Vaucouleurs (1948) (n=4) profiles.
 *
 * This particular implementation closely matches the majority of the profile,
 * but deliberately excludes the very inner and outermost regions.
 *
 * @note See [Graham & Driver 2005](https://ned.ipac.caltech.edu/level5/March05/Graham/Graham2.html)
 *       for a useful summary of various properties of the Sersic profile.
 */
class SersicMixComponent : private SersicParametricEllipseHolder, public EllipticalComponent {
public:
    explicit SersicMixComponent(std::shared_ptr<SersicParametricEllipse> ellipse = nullptr,
                                std::shared_ptr<CentroidParameters> centroid = nullptr,
                                std::shared_ptr<IntegralModel> integralmodel = nullptr,
                                std::shared_ptr<SersicMixComponentIndexParameterD> sersicindex = nullptr);
    ~SersicMixComponent();

    void add_extra_param_map(const Channel& channel, ExtraParamMap& map_extra, const GradParamMap& map_grad,
                             ParameterMap& offsets) const override;
    void add_extra_param_factors(const Channel& channel, ExtraParamFactors& factors) const override;
    void add_grad_param_map(const Channel& channel, GradParamMap& map, ParameterMap& offsets) const override;
    void add_grad_param_factors(const Channel& channel, GradParamFactors& factor) const override;

    std::unique_ptr<const lsst::gauss2d::Gaussians> get_gaussians(const Channel& channel) const override;
    size_t get_n_gaussians(const Channel& channel) const override;

    ParamRefs& get_parameters(ParamRefs& params, ParamFilter* filter = nullptr) const override;
    ParamCRefs& get_parameters_const(ParamCRefs& params, ParamFilter* filter = nullptr) const override;

    double get_sersicindex() const;

    SersicMixComponentIndexParameterD& get_sersicindex_param() const;
    std::shared_ptr<SersicMixComponentIndexParameterD> get_sersicindex_param_ptr();

    static const size_t N_PARAMS = N_PARAMS_GAUSS2D + 1;

    void set_extra_param_factors(const Channel& channel, ExtraParamFactors& factors,
                                 size_t index) const override;
    void set_grad_param_factors(const Channel& channel, GradParamFactors& factors,
                                size_t index) const override;
    void set_sersicindex(double value);

    std::string repr(bool name_keywords = false,
                     std::string_view namespace_separator = Object::CC_NAMESPACE_SEPARATOR) const override;
    std::string str() const override;

private:
    class SersicMixGaussianComponent;
    std::shared_ptr<SersicMixComponentIndexParameterD> _sersicindex;
    std::map<std::reference_wrapper<const Channel>, std::vector<std::unique_ptr<SersicMixGaussianComponent>>>
            _gaussians;
};
}  // namespace lsst::gauss2d::fit

#endif
