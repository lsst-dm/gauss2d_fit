#include <stdexcept>

#include <optional>
#include "source.h"

namespace gauss2d
{
namespace fit
{
void Source::add_extra_param_map(
    const Channel & channel, extra_param_map & map_extra,
    const grad_param_map & map_grad, ParameterMap & offsets
) const
{
    for(auto & component : _components) component->add_extra_param_map(
        channel, map_extra, map_grad, offsets);
}

void Source::add_extra_param_factors(const Channel & channel, extra_param_factors & factors) const
{
    for(auto & component : _components) component->add_extra_param_factors(channel, factors);
}

void Source::add_grad_param_map(const Channel & channel, grad_param_map & map, ParameterMap & offsets
    ) const
{
    for(auto & component : _components) component->add_grad_param_map(channel, map, offsets);
}

void Source::add_grad_param_factors(const Channel & channel, grad_param_factors & factors) const
{
    for(auto & component : _components) component->add_grad_param_factors(channel, factors);
}

std::unique_ptr<const gauss2d::Gaussians> Source::get_gaussians(const Channel & channel) const 
{
    std::vector<std::optional<const gauss2d::Gaussians::Data>> in;
    // TODO: This isn't sufficient; need to implement get_n_components
    in.reserve(_components.size());
    for(auto & component : _components) {
        in.push_back(component->get_gaussians(channel)->get_data());
    }
    return std::make_unique<gauss2d::Gaussians>(in);
}

size_t Source::get_n_gaussians(const Channel & channel) const {
    size_t n_g = 0;
    for(auto & component : _components) n_g += component->get_n_gaussians(channel);
    return n_g;
}

ParamRefs & Source::get_parameters(ParamRefs & params, ParamFilter * filter) const
{
    for(auto & component : _components) component->get_parameters(params, filter);
    return params;
}

ParamCRefs & Source::get_parameters_const(ParamCRefs & params, ParamFilter * filter) const
{
    for(auto & component : _components) component->get_parameters_const(params, filter);
    return params;
}

void Source::set_extra_param_factors(const Channel & channel, extra_param_factors & factors, size_t index) const {
    for(auto & component : _components) {
        component->set_extra_param_factors(channel, factors, index);
        index += component->get_n_gaussians(channel);
    }
}

void Source::set_grad_param_factors(const Channel & channel, grad_param_factors & factors, size_t index) const {
    for(auto & component : _components) {
        component->set_grad_param_factors(channel, factors, index);
        index += component->get_n_gaussians(channel);
    }
}

std::string Source::str() const {
    std::string s = "Source([";
    for(const auto & c : _components) s += c->str() + ",";
    return s + "])";
}

Source::Source(Source::Components & components)
{
    _components.reserve(components.size());
    size_t i = 0;
    for(auto & component : components) {
        if(component == nullptr) throw std::invalid_argument("Source components[" + std::to_string(i)
            + "] can't be null");
        _components.push_back(std::move(component));
        i++;
    }
}

} // namespace fit
} // namespace gauss2d
