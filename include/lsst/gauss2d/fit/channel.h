#ifndef LSST_GAUSS2D_FIT_CHANNEL_H
#define LSST_GAUSS2D_FIT_CHANNEL_H

#include <iostream>
#include <functional>
#include <map>
#include <memory>
#include <set>
#include <stdexcept>
#include <string>

#include "lsst/gauss2d/object.h"

#include "util.h"

namespace lsst::gauss2d::fit {

/**
 * @brief An observational channel, usually representing some range of wavelengths of light.
 *
 * A Channel is a distinguishing property of an Observation, allowing
 * IntegralModels to define integrals in distinct regions of some parameter
 * space. Channels can represent a physical object such as a filter or more
 * abstract definitions like wavelength/frequency ranges, generic labels
 * like ABC, something inbetween like RGB, or anything else users prefer.
 *
 * @note Channels are compared and sorted alphabetically.
 */
class Channel : public Object {
public:
    typedef std::map<std::string, std::shared_ptr<const Channel>> Registry;

    Channel(const Channel &) = delete;
    Channel &operator=(const Channel &) = delete;
    /**
     * Delete a channel with a given name.
     *
     * @param name The name of the channel to delete
     *
     * @throw std::invalid_argument If the channel is the default (NAME_NONE) or non-existent
     * @throw std::runtime_error If the channel is referenced by any other objects
     *
     */
    static void erase(std::string name);

    /**
     * Find a channel with the given name
     *
     * @param name The name of the channel to find
     * @return The Channel is if it exists, or nullptr otherwise
     */
    static const std::shared_ptr<const Channel> find_channel(std::string name);

    /**
     * Get a channel with the given name, creating it if it doesn't exist
     *
     * @param name The name of the channel to get
     * @return The Channel with the given name (created if non-existent)
     */
    static const std::shared_ptr<const Channel> get_channel(std::string name);

    static std::vector<std::shared_ptr<const Channel>> get_channels();

    inline static const std::string NAME_NONE = "None";

    const std::string name;

    static const std::shared_ptr<const Channel> NONE_PTR();
    static const Channel &NONE();

    std::string repr(bool name_keywords = false,
                     std::string_view namespace_separator = Object::CC_NAMESPACE_SEPARATOR) const override;
    std::string str() const override;

    /**
     * Construct a new Channel
     *
     * @param name The name of the Channel to create
     * @return The new Channel
     */
    static std::shared_ptr<Channel> make(std::string name);
    /// Same as make(), but creating a new Channel
    static const std::shared_ptr<const Channel> make_const(std::string name);

    // TODO: Figure out why tests compile but do not run without this operator.
    // Until then, do NOT remove.
    const bool operator<(const Channel &c) const;
    const bool operator==(const Channel &c) const;
    const bool operator!=(const Channel &c) const;

private:
    /**
     * @brief Channel's private constructor.
     *
     * This constructor is private to force all newly-made Channels to be registered.
     *
     * @param name The name of the channel
     */
    explicit Channel(std::string name);

    struct Shared_enabler;
};

inline bool operator<(const std::reference_wrapper<const Channel> &lhs,
                      const std::reference_wrapper<const Channel> &rhs) {
    return (lhs.get() < rhs.get());
}

inline bool operator==(const std::reference_wrapper<const Channel> &lhs,
                       const std::reference_wrapper<const Channel> &rhs) {
    return (lhs.get() == rhs.get());
}

}  // namespace lsst::gauss2d::fit

#endif
