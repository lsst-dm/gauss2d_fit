#include "psfmodel.h"

#include <memory>
#include <optional>

#include "param_filter.h"
#include "source.h"

namespace gauss2d
{
namespace fit
{
void PsfModel::add_extra_param_map(
    const Channel & channel, extra_param_map & map_extra,
    const grad_param_map & map_grad, ParameterMap & offsets
) const
{
    for(auto & component : _components) component->add_extra_param_map(
        channel, map_extra, map_grad, offsets);
}

void PsfModel::add_extra_param_factors(const Channel & channel, extra_param_factors & factors) const
{
    for(auto & component : _components) component->add_extra_param_factors(channel, factors);
}

void PsfModel::add_grad_param_map(const Channel & channel, grad_param_map & map, ParameterMap & offsets
    ) const
{
    for(auto & component : _components) component->add_grad_param_map(channel, map, offsets);
}

void PsfModel::add_grad_param_factors(const Channel & channel, grad_param_factors & factors) const
{
    for(auto & component : _components) component->add_grad_param_factors(channel, factors);
}

Components PsfModel::get_components() const {
    return _components;
}

std::unique_ptr<const gauss2d::Gaussians> PsfModel::get_gaussians(const Channel & channel) const 
{
    std::vector<std::optional<const gauss2d::Gaussians::Data>> in;
    in.reserve(_components.size());
    for(auto & component : _components) {
        in.push_back(component->get_gaussians(channel)->get_data());
    }
    return std::make_unique<gauss2d::Gaussians>(in);
}

size_t PsfModel::get_n_gaussians(const Channel & channel) const {
    size_t n_g = 0;
    for(auto & component : _components) n_g += component->get_n_gaussians(channel);
    return n_g;
}

ParamRefs & PsfModel::get_parameters(ParamRefs & params, ParamFilter * filter) const
{
    for(auto & component : _components) component->get_parameters(params, filter);
    return params;
}
ParamCRefs & PsfModel::get_parameters_const(ParamCRefs & params, ParamFilter * filter) const
{
    for(auto & component : _components) component->get_parameters_const(params, filter);
    return params;
}

void PsfModel::set_extra_param_factors(const Channel & channel, extra_param_factors & factors, size_t index) const
{
    for(auto & component : _components) component->set_extra_param_factors(channel, factors, index);
}

void PsfModel::set_grad_param_factors(const Channel & channel, grad_param_factors & factors, size_t index) const
{
    for(auto & component : _components) component->set_grad_param_factors(channel, factors, index);
}

std::string PsfModel::str() const {
    std::string str = "PsfModel(components=[";
    for(const auto & s : _components) str += s->str() + ",";
    return str + "])";
}

PsfModel::PsfModel(Components & components) {
    size_t i = 0;
    _components.reserve(components.size());
    for(auto & component : components) {
       if(component == nullptr) throw std::invalid_argument("PsfModel components["
            + std::to_string(i) + "] can't be null");
        auto channels = component->get_integralmodel().get_channels();
        if((channels.size() != 1) || ((*channels.begin()).get() != Channel::NONE())) {
            throw std::invalid_argument("PsfModel components[" + std::to_string(i)
                + "].get_integralmodel().get_channels()=" + str_iter_refw(channels)
                + " must only contain None");
        }
        _components.push_back(std::move(component));
        i++;
    }
}

PsfModel::~PsfModel() {}
} // namespace fit
} // namespace gauss2d
