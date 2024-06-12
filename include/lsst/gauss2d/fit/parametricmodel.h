#ifndef LSST_GAUSS2D_FIT_PARAMETRICMODEL_H
#define LSST_GAUSS2D_FIT_PARAMETRICMODEL_H

#include "lsst/gauss2d/evaluate.h"
#include "lsst/gauss2d/gaussian.h"

#include "channel.h"
#include "parametric.h"

namespace lsst::gauss2d::fit{
typedef std::vector<std::array<size_t, lsst::gauss2d::N_EXTRA_MAP>> ExtraParamMap;
typedef std::array<double, lsst::gauss2d::N_EXTRA_FACTOR> ExtraParamFactorValues;
typedef std::vector<ExtraParamFactorValues> ExtraParamFactors;
typedef std::vector<std::array<size_t, lsst::gauss2d::N_PARAMS_GAUSS2D>> GradParamMap;
typedef std::vector<std::array<double, lsst::gauss2d::N_PARAMS_GAUSS2D>> GradParamFactors;

typedef std::map<ParamBaseCRef, size_t> ParameterMap;

// TODO: Expand on the method docs; add reasoning on the design
/**
 * A Parametric that can manage Parameter gradients and return Gaussians.
 *
 * @note See GaussianEvaluator for details on extra/grad param maps.
 */
class ParametricModel : public Parametric {
public:
    /**
     * Add extra Parameter indices to a map.
     *
     * @param channel The Channel to add indices for.
     * @param map_extra The ExtraParamMap to add to.
     * @param map_grad The completed GradParamMap.
     * @param offsets A map of index offsets for Parameters that have already been added.
     */
    virtual void add_extra_param_map(const Channel& channel, ExtraParamMap& map_extra,
                                     const GradParamMap& map_grad, ParameterMap& offsets) const = 0;
    /**
     * Add extra Parameter gradient factors to an existing vector.
     *
     * @param channel The Channel to add factors for.
     * @param factors The ExtraParamFactors to add to.
     */
    virtual void add_extra_param_factors(const Channel& channel, ExtraParamFactors& factors) const = 0;
    /**
     * Add Parameter gradient indices to an existing map.
     *
     * @param channel The Channel to add indices for.
     * @param map The map to add to.
     * @param offsets A map of index offsets for Parameters that have already been added.
     */
    virtual void add_grad_param_map(const Channel& channel, GradParamMap& map,
                                    ParameterMap& offsets) const = 0;
    /**
     * Add Parameter gradient factors to an existing map.
     *
     * @param channel The Channel to add factors for.
     * @param factors The GradParamFactors to add to.
     */
    virtual void add_grad_param_factors(const Channel& channel, GradParamFactors& factors) const = 0;
    /**
     * Set extra Parameter gradient factors in an existing map.
     *
     * @param channel The Channel to set factors for.
     * @param factors The ExtraParamFactors to set factors for.
     * @param index The index to begin setting factors at.
     */
    virtual void set_extra_param_factors(const Channel& channel, ExtraParamFactors& factors,
                                         size_t index) const = 0;
    /**
     * Set Parameter gradient factors in an existing map.
     *
     * @param channel The Channel to set factors for.
     * @param factors The GradParamFactors to set factors for.
     * @param index The index to begin setting factors at.
     */
    virtual void set_grad_param_factors(const Channel& channel, GradParamFactors& factors,
                                        size_t index) const = 0;

    /**
     * Return the vector of Gaussian sub-components controlled by this model.
     *
     * @param channel The Channel to return Gaussians for.
     * @return The Gaussians controlled by this model.
     */
    virtual std::unique_ptr<const lsst::gauss2d::Gaussians> get_gaussians(const Channel& channel) const = 0;
    /// Return the number of Gaussian sub-components controlled by this model.
    virtual size_t get_n_gaussians(const Channel& channel) const = 0;
};

}  // namespace lsst::gauss2d::fit

#endif
