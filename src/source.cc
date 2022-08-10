#include <stdexcept>

#include <optional>
#include "source.h"

namespace gauss2d
{
namespace fit
{
std::unique_ptr<const gauss2d::Gaussians> Source::get_gaussians(const Channel & channel) const 
{
    std::vector<std::optional<const gauss2d::Gaussians::Data>> in;
    in.reserve(_components.size());
    for(auto & component : _components) {
        in.push_back(component->get_gaussians(channel)->get_data());
    }
    return std::make_unique<gauss2d::Gaussians>(in);
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
