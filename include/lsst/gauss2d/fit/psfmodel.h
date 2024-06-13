#ifndef LSST_GAUSS2D_FIT_PSFMODEL_H
#define LSST_GAUSS2D_FIT_PSFMODEL_H

#include <memory>

#include "component.h"
#include "componentmixture.h"
#include "param_filter.h"

namespace lsst::gauss2d::fit {
/**
 * @brief A Gaussian mixture model of a point spread function.
 *
 * A PsfModel is, like a Source, a ComponentMixture. It represents the PSF,
 * i.e. the smoothing kernel for a single Observation (whether from the
 * optical system, environmental conditions, or any other source of blurring).
 * As such, it should have an IntegralModel instance that sum to unity.
 * This is most easily enforced with the use of FractionalIntegralModel.
 *
 * PsfModels are also generally required to not have a specific Channel.
 * Logically, it should have the same Channel as the Observation it applies
 * to, but generally, it cannot be defined to apply to multiple Observations,
 * so non-NONE Channels are disallowed to reflect this fact.
 *
 */
class PsfModel : public ComponentMixture {
public:
    explicit PsfModel(Components& components);
    ~PsfModel();
    
    void add_extra_param_map(const Channel& channel, ExtraParamMap& map_extra, const GradParamMap& map_grad,
                             ParameterMap& offsets) const override;
    void add_extra_param_factors(const Channel& channel, ExtraParamFactors& factors) const override;
    void add_grad_param_map(const Channel& channel, GradParamMap& map, ParameterMap& offsets) const override;
    void add_grad_param_factors(const Channel& channel, GradParamFactors& factor) const override;

    Components get_components() const override;
    std::unique_ptr<const lsst::gauss2d::Gaussians> get_gaussians(const Channel& channel
                                                                  = Channel::NONE()) const override;
    size_t get_n_gaussians(const Channel& channel = Channel::NONE()) const override;

    ParamRefs& get_parameters(ParamRefs& params, ParamFilter* filter = nullptr) const override;
    ParamCRefs& get_parameters_const(ParamCRefs& params, ParamFilter* filter = nullptr) const override;

    void set_extra_param_factors(const Channel& channel, ExtraParamFactors& factors,
                                 size_t index) const override;
    void set_grad_param_factors(const Channel& channel, GradParamFactors& factor,
                                size_t index) const override;

    std::string repr(bool name_keywords = false,
                     std::string_view namespace_separator = Object::CC_NAMESPACE_SEPARATOR) const override;
    std::string str() const override;

private:
    Components _components = {};
};

}  // namespace lsst::gauss2d::fit

#endif
