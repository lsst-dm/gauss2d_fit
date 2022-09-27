#include "psfmodel.h"

#include <memory>
#include <optional>

#include "param_filter.h"
#include "source.h"

namespace gauss2d
{
namespace fit
{
std::unique_ptr<const gauss2d::Gaussians> PsfModel::get_gaussians(const Channel & channel) const 
{
    std::vector<std::optional<const gauss2d::Gaussians::Data>> in;
    in.reserve(_components.size());
    for(auto & component : _components) {
        in.push_back(component->get_gaussians(channel)->get_data());
    }
    return std::make_unique<gauss2d::Gaussians>(in);
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
