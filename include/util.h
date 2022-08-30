#ifndef GAUSS2D_FIT_UTIL_H
#define GAUSS2D_FIT_UTIL_H

#include <memory>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>

namespace gauss2d
{
namespace fit
{

// following numpy.isclose
template <typename T>
struct IsCloseResult
{
    const bool isclose;
    const double diff_abs;
    const T margin;

    std::string str() const {
        std::stringstream ss;
        ss.precision(5);
        ss << "isclose=" << isclose << " from diff=" << diff_abs << " <= margin=" << margin;
        return ss.str();
    }
};

template <typename T>
IsCloseResult<T> isclose(T a, T b, T rtol=1e-5, T atol=1e-8)
{
    const T diff_abs = std::abs(a - b);
    const T margin = atol + rtol*std::abs(b);
    return IsCloseResult<T>{diff_abs <= margin, diff_abs, margin};
}

// The rest of these functions are mainly intended to make printing container members easier
template<template <typename...> class Map, class Key, class Value>
std::set<Key>
map_keys(const Map<Key, Value> & map)
{
    std::set<Key> result;
    for(const auto& it : map){
        result.insert(it.first);
    }
    return result;
}

template<template <typename...> class Map, class Key, class Value>
std::set<std::reference_wrapper<const Key>>
map_keys_ref_const(const Map<std::reference_wrapper<const Key>, Value, std::less<const Key>> & map)
{
    std::set<std::reference_wrapper<const Key>> result;
    for(const auto& it : map) result.insert(it.first.get());
    return result;
}

template<template <typename...> class Map, class Key, class Value>
std::set<std::shared_ptr<const Value>>
map_values_shared_ptr_const(
    const Map<Key, std::shared_ptr<const Value>> & map
)
{
    std::set<std::shared_ptr<const Value>> result;
    for(const auto& it : map) result.insert(it.second);
    return result;
}

template<template <typename...> class Map, class Key, class Value>
std::set<std::reference_wrapper<const Value>>
map_values_ref_const(
    const Map<Key, std::reference_wrapper<const Value>> & map
)
{
    std::set<std::reference_wrapper<const Value>> result;
    for(const auto& it : map) result.insert(it.second);
    return result;
}

template <typename T>
std::string str_ptr(const T* obj)
{
    return (obj == nullptr) ? "None" : obj->str();
}

template <typename T>
std::string str_iter_ptr(const T & container)
{
    std::string str = "[";
    for(const auto & obj: container) str += str_ptr(obj) + ",";
    return str.substr(0, str.size() - 1) + "]";
}

template <typename T>
std::string str_iter_ref(const T & container)
{
    std::string str = "[";
    for(const auto & obj: container) str += obj.str() + ",";
    return str.substr(0, str.size() - 1) + "]";
}

template <typename T>
std::string str_iter_refw(const T & container)
{
    std::string str = "[";
    for(const auto & obj: container) str += obj.get().str() + ",";
    return str.substr(0, str.size() - 1) + "]";
}

} // namespace fit
} // namespace gauss2d

#endif
