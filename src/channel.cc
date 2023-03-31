#include "channel.h"

#include <functional>
#include <string>

#include "util.h"

#include <iostream>

namespace gauss2d
{
namespace fit
{

static inline Channel::Registry _registry = {};

// https://stackoverflow.com/questions/8147027/
// how-do-i-call-stdmake-shared-on-a-class-with-only-protected-or-private-const/
// 8147213#comment58654091_25069711
struct Channel::Shared_enabler : public Channel
{
    template <typename... Args>
    Shared_enabler(Args &&... args)
    : Channel(std::forward<Args>(args)...) {}
};

const bool Channel::operator < ( const Channel &c ) const {
    return name < c.name;
}

const bool Channel::operator == ( const Channel &c ) const {
    return name == c.name;
}

const bool Channel::operator != ( const Channel &c ) const {
    return !(*this == c);
}

void Channel::erase(std::string name) {
    if(name == NAME_NONE) throw std::invalid_argument("Can't erase the " + NAME_NONE + " Channel");
    const std::shared_ptr<const Channel> channel = Channel::get_channel(name);
    if(channel == nullptr) throw std::invalid_argument("Can't erase non-existent Channel name=" + name);
    size_t use_count = channel.use_count();
    if(use_count < 2) {
        throw std::logic_error("Failed erasing " + channel->str() + " with unexpected use_count="
            +  std::to_string(use_count));
    } else if(use_count == 2) {
        _registry.erase(name);
    } else {
        throw std::runtime_error("Can't erase " + channel->str() + " with use_count="
            + std::to_string(use_count));
    }
}

const std::shared_ptr<const Channel> Channel::find_channel(std::string name) {
    if(name == NAME_NONE) return NONE_PTR();
    auto channel = _registry.find(name);
    auto found = channel == _registry.end() ? nullptr : channel->second;
    return found;
}

const std::shared_ptr<const Channel> Channel::get_channel(std::string name) {
    auto channel = Channel::find_channel(name);
    if(channel == nullptr) channel = Channel::make(name);
    return channel;
}

std::set<std::shared_ptr<const Channel>> Channel::get_channels() {
    auto set = map_values_shared_ptr_const(_registry);
    set.insert(NONE_PTR());
    return set;
}

std::shared_ptr<Channel> Channel::make(std::string name) {
    std::shared_ptr<Channel> c = std::make_shared<Shared_enabler>(name);
    if(name != NAME_NONE) {
        _registry.insert({name, c});
    }
    return c;
}

const std::shared_ptr<const Channel> Channel::make_const(std::string name) {
    return make(name);
}

std::string Channel::repr(bool name_keywords) const {
    return std::string("Channel(") + (name_keywords ? "name=" : "") + name + ")";
}

std::string Channel::str() const { return "Channel(name=" + name + ")"; }

const std::shared_ptr<const Channel> Channel::NONE_PTR() {
    static const std::shared_ptr<const Channel> _NONE = Channel::make(Channel::NAME_NONE);
    return _NONE;
}

const Channel & Channel::NONE() {
    return *NONE_PTR();
}

Channel::Channel(std::string name_) : name(name_) {
    static bool none_made = false;
    if(name_ == NAME_NONE) {
        if(none_made) {
            throw std::invalid_argument("Channel('None') always exists and cannot be re-made");
        } else {
            none_made = true;
        }
    } else {
        auto found = _registry.find(name);
        if(found != _registry.end()) {
            throw std::invalid_argument("Channel(name=" + name + ") already exists as "
                + found->second->str()
            );
        }
    }
}

} // namespace fit
} // namespace gauss2d
