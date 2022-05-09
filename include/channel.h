#ifndef GAUSS2D_FIT_CHANNEL_H
#define GAUSS2D_FIT_CHANNEL_H

#include <iostream>
#include <functional>
#include <map>
#include <memory>
#include <set>
#include <stdexcept>
#include <string>

#include "gauss2d/object.h"

#include "util.h"

namespace gauss2d
{
namespace fit
{

class Channel : public Object//, public std::enable_shared_from_this<Channel>
{
public:
    typedef std::map<std::string, std::shared_ptr<const Channel>> Registry;

private:
    static inline Registry _registry = {};
    Channel(std::string name);

    struct Shared_enabler;

public:
    static void erase(std::string name);

    static const std::shared_ptr<const Channel> get_channel(std::string name) {
        if(name == NAME_NONE) return NONE_PTR();
        auto channel = _registry.find(name);
        return channel == _registry.end() ? nullptr : channel->second;
    }
    static std::set<std::shared_ptr<const Channel>> get_channels() {
        auto set = map_values_shared_ptr_const(_registry);
        set.insert(NONE_PTR());
        return set;
    }

    inline static const std::string NAME_NONE = "None";

    const std::string name;

    static const std::shared_ptr<const Channel> NONE_PTR();
    static const Channel & NONE();

    std::string str() const override;

    static std::shared_ptr<Channel> make(std::string name);
    static const std::shared_ptr<const Channel> make_const(std::string name);

    // TODO: Figure out why tests compile but do not run without this operator.
    // Until then, do NOT remove.
    const bool operator < ( const Channel &c ) const;
    const bool operator == ( const Channel &c ) const;
    const bool operator != ( const Channel &c ) const;

    Channel (const Channel&) = delete;
    Channel& operator= (const Channel&) = delete;
};

inline bool operator < (const std::reference_wrapper<const Channel> & lhs,
                 const std::reference_wrapper<const Channel> & rhs) {
    return (lhs.get() < rhs.get());
}

inline bool operator == (const std::reference_wrapper<const Channel> & lhs,
                 const std::reference_wrapper<const Channel> & rhs) {
    return (lhs.get() == rhs.get());
}

} // namespace fit
} // namespace gauss2d

#endif