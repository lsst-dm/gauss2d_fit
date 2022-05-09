#ifndef GAUSS2D_FIT_UTIL_H
#define GAUSS2D_FIT_UTIL_H

#include <memory>
#include <set>
#include <stdexcept>
#include <string>

namespace gauss2d
{
namespace fit
{

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