#ifndef LSST_GAUSS2D_FIT_GAUSSIANCOMPONENT_H
#define LSST_GAUSS2D_FIT_GAUSSIANCOMPONENT_H

#include "channel.h"
#include "ellipticalcomponent.h"
#include "gaussianparametricellipse.h"
#include "integralmodel.h"
#include "linearintegralmodel.h"
#include "param_defs.h"
#include "param_filter.h"
#include <memory>

namespace lsst::gauss2d::fit {
// TODO: Revisit the necessity of this class
// Its purpose is to have the GaussianParametricEllipse stored here and initialized first in
// GaussianComponent's constructor
class GaussianParametricEllipseHolder {
public:
    std::shared_ptr<GaussianParametricEllipse> _ellipsedata;

    explicit GaussianParametricEllipseHolder(std::shared_ptr<GaussianParametricEllipse> ellipse = nullptr)
            : _ellipsedata(std::move(ellipse)) {
        if (_ellipsedata == nullptr) _ellipsedata = std::make_shared<GaussianParametricEllipse>();
    }
};

/**
 * A Component consisting of a 2D Gaussian.
 */
class GaussianComponent : private GaussianParametricEllipseHolder, public EllipticalComponent {
public:
    /**
     * Construct a GaussianComponent from ellipse, centroid and integral parameters.
     *
     * @param ellipse The GaussianParametricEllipse value; default-initialized if null.
     * @param centroid The CentroidParameters value; default-initialized if null.
     * @param integralmodel The IntegralModel value; default-initialized if null.
     */
    explicit GaussianComponent(std::shared_ptr<GaussianParametricEllipse> ellipse = nullptr,
                               std::shared_ptr<CentroidParameters> centroid = nullptr,
                               std::shared_ptr<IntegralModel> integralmodel = nullptr);

    void add_extra_param_map(const Channel& channel, ExtraParamMap& map_extra, const GradParamMap& map_grad,
                             ParameterMap& offsets) const override;
    void add_extra_param_factors(const Channel& channel, ExtraParamFactors& factors) const override;
    void add_grad_param_map(const Channel& channel, GradParamMap& map, ParameterMap& offsets) const override;
    void add_grad_param_factors(const Channel& channel, GradParamFactors& factor) const override;

    std::unique_ptr<const lsst::gauss2d::Gaussians> get_gaussians(const Channel& channel) const override;
    size_t get_n_gaussians(const Channel& channel) const override;

    ParamRefs& get_parameters(ParamRefs& params, ParamFilter* filter = nullptr) const override;
    ParamCRefs& get_parameters_const(ParamCRefs& params, ParamFilter* filter = nullptr) const override;

    /**
     * Construct a vector of default-initialized GaussianComponent instances.
     *
     * @param sizes Vector of initial values for both sigma_x and sigma_y.
     * @param fixed Whether all Parameter members should be fixed initially.
     * @return A vector of GaussianComponent instances.
     *
     * @note This can be used to initialize a trivial PsfModel as a
     * noralized, zero-size single Gaussian.
     */
    static std::vector<std::shared_ptr<Component>> make_uniq_default_gaussians(
            const std::vector<double>& sizes = {2.}, bool fixed = true) {
        std::vector<std::shared_ptr<Component>> comps = {};
        for (const double size : sizes) {
            LinearIntegralModel::Data data
                    = {{Channel::NONE(),
                        std::make_shared<IntegralParameterD>(1., nullptr, nullptr, nullptr, fixed)}};
            comps.emplace_back(std::make_shared<GaussianComponent>(
                    std::make_shared<GaussianParametricEllipse>(
                            std::make_shared<SigmaXParameterD>(size, nullptr, nullptr, nullptr, fixed),
                            std::make_shared<SigmaYParameterD>(size, nullptr, nullptr, nullptr, fixed),
                            std::make_shared<RhoParameterD>(0, nullptr, nullptr, nullptr, fixed)),
                    std::make_shared<CentroidParameters>(
                            std::make_shared<CentroidXParameterD>(0, nullptr, nullptr, nullptr, fixed),
                            std::make_shared<CentroidYParameterD>(0, nullptr, nullptr, nullptr, fixed)),
                    std::make_shared<LinearIntegralModel>(&data)));
        }
        return comps;
    }

    void set_extra_param_factors(const Channel& channel, ExtraParamFactors& factors,
                                 size_t index) const override;
    void set_grad_param_factors(const Channel& channel, GradParamFactors& factors,
                                size_t index) const override;

    std::string repr(bool name_keywords = false,
                     std::string_view namespace_separator = Object::CC_NAMESPACE_SEPARATOR) const override;
    std::string str() const override;

private:
    ParamCRefs _get_parameters_grad(const Channel& channel) const;
};
}  // namespace lsst::gauss2d::fit

#endif
