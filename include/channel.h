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

/*
    A Channel is a distinguishing property of an Observation, allowing
    IntegralModels to define integrals in distinct regions of some parameter
    space. Channels can represent a physical object such as a filter or more
    abstract definitions like wavelength/frequency ranges, generic labels
    like ABC, something inbetween like RGB, or anything else users prefer.
*/
class Channel : public Object
{
public:
    typedef std::map<std::string, std::shared_ptr<const Channel>> Registry;

private:
/*
    The registry serves several purposes:
        - ensures that Channels are unique and do not share names
        - allows Channels to be easily found by name.
        - prevents Channels from being implicitly deleted when no objects reference them
          (though they can be explicitly deleted in that case)
*/
    Channel(std::string name);

    struct Shared_enabler;

public:
    static void erase(std::string name);

    static const std::shared_ptr<const Channel> find_channel(std::string name);

    static const std::shared_ptr<const Channel> get_channel(std::string name);

    static std::set<std::shared_ptr<const Channel>> get_channels();

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
