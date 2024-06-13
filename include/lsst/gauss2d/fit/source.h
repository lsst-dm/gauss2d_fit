#ifndef LSST_GAUSS2D_FIT_SOURCE_H
#define LSST_GAUSS2D_FIT_SOURCE_H

#include <memory>

#include "component.h"
#include "componentmixture.h"
#include "param_filter.h"

namespace lsst::gauss2d::fit {
/**
 * @brief An association of physically-related Component instances.
 *
 * A Source is primarily intended to be a single object represented by
 * Component instances with common centroids. However, users may link
 * whatever Components they like into a Source.
 *
 * Source objects should not share Components if they are part of the
 * same Model, but this is not enforced.
 *
 */
class Source : public ComponentMixture {
public:
    explicit Source(Components& components);

    void add_extra_param_map(const Channel& channel, ExtraParamMap& map_extra, const GradParamMap& map_grad,
                             ParameterMap& offsets) const override;
    void add_extra_param_factors(const Channel& channel, ExtraParamFactors& factors) const override;
    void add_grad_param_map(const Channel& channel, GradParamMap& map, ParameterMap& offsets) const override;
    void add_grad_param_factors(const Channel& channel, GradParamFactors& factors) const override;

    Components get_components() const override;
    std::unique_ptr<const lsst::gauss2d::Gaussians> get_gaussians(const Channel& channel) const override;
    size_t get_n_gaussians(const Channel& channel) const override;

    ParamRefs& get_parameters(ParamRefs& params, ParamFilter* filter = nullptr) const override;
    ParamCRefs& get_parameters_const(ParamCRefs& params, ParamFilter* filter = nullptr) const override;

    void set_extra_param_factors(const Channel& channel, ExtraParamFactors& factors,
                                 size_t index) const override;
    void set_grad_param_factors(const Channel& channel, GradParamFactors& factors,
                                size_t index) const override;

    std::string repr(bool name_keywords = false,
                     std::string_view namespace_separator = Object::CC_NAMESPACE_SEPARATOR) const override;
    std::string str() const override;

private:
    Components _components = {};
};

}  // namespace lsst::gauss2d::fit

#endif
