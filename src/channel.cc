#include <functional>
#include <iostream>
#include <string>

#include "lsst/gauss2d/type_name.h"

#include "lsst/gauss2d/fit/channel.h"

namespace lsst::gauss2d::fit {

/*
The registry serves several purposes:
    - ensures that Channels are unique and do not share names
    - allows Channels to be easily found by name.
    - prevents Channels from being implicitly deleted when no objects reference them
      (though they can be explicitly deleted in that case)
*/
static inline Channel::Registry _registry = {};

// https://stackoverflow.com/questions/8147027/
// how-do-i-call-stdmake-shared-on-a-class-with-only-protected-or-private-const/
// 8147213#comment58654091_25069711
struct Channel::Shared_enabler : public Channel {
    template <typename... Args>
    Shared_enabler(Args &&...args) : Channel(std::forward<Args>(args)...) {}
};

Channel::Channel(std::string name_) : name(name_) {
    static bool none_made = false;
    if (name_ == NAME_NONE) {
        if (none_made) {
            throw std::invalid_argument("Channel('None') always exists and cannot be re-made");
        } else {
            none_made = true;
        }
    } else {
        auto found = _registry.find(name);
        if (found != _registry.end()) {
            throw std::invalid_argument("Channel(name=" + name + ") already exists as "
                                        + found->second->str());
        }
    }
}

const bool Channel::operator<(const Channel &c) const { return name < c.name; }

const bool Channel::operator==(const Channel &c) const { return name == c.name; }

const bool Channel::operator!=(const Channel &c) const { return !(*this == c); }

void Channel::erase(std::string name) {
    if (name == NAME_NONE) throw std::invalid_argument("Can't erase the " + NAME_NONE + " Channel");
    const std::shared_ptr<const Channel> channel = Channel::get_channel(name);
    if (channel == nullptr) throw std::invalid_argument("Can't erase non-existent Channel name=" + name);
    size_t use_count = channel.use_count();
    if (use_count < 2) {
        throw std::logic_error("Failed erasing " + channel->str()
                               + " with unexpected use_count=" + std::to_string(use_count));
    } else if (use_count == 2) {
        _registry.erase(name);
    } else {
        throw std::runtime_error("Can't erase " + channel->str()
                                 + " with use_count=" + std::to_string(use_count));
    }
}

const std::shared_ptr<const Channel> Channel::find_channel(std::string name) {
    if (name == NAME_NONE) return NONE_PTR();
    auto channel = _registry.find(name);
    auto found = channel == _registry.end() ? nullptr : channel->second;
    return found;
}

const std::shared_ptr<const Channel> Channel::get_channel(std::string name) {
    auto channel = Channel::find_channel(name);
    if (channel == nullptr) channel = Channel::make(name);
    return channel;
}

std::vector<std::shared_ptr<const Channel>> Channel::get_channels() {
    std::vector<std::shared_ptr<const Channel>> vec = {};
    for (const auto &[key, value] : _registry) {
        vec.emplace_back(value);
    }
    vec.emplace_back(NONE_PTR());
    return vec;
}

std::shared_ptr<Channel> Channel::make(std::string name) {
    std::shared_ptr<Channel> c = std::make_shared<Shared_enabler>(name);
    if (name != NAME_NONE) {
        _registry.insert({name, c});
    }
    return c;
}

const std::shared_ptr<const Channel> Channel::make_const(std::string name) { return make(name); }

std::string Channel::repr(bool name_keywords, std::string_view namespace_separator) const {
    return type_name_str<Channel>(false, namespace_separator) + "(" + (name_keywords ? "name=" : "") + name
           + ")";
}

std::string Channel::str() const { return type_name_str<Channel>(true) + "(name=" + name + ")"; }

const std::shared_ptr<const Channel> Channel::NONE_PTR() {
    static const std::shared_ptr<const Channel> _NONE = Channel::make(Channel::NAME_NONE);
    return _NONE;
}

const Channel &Channel::NONE() { return *NONE_PTR(); }

}  // namespace lsst::gauss2d::fit
