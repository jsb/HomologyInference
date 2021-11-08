/*
 * Author: Patrick Schmidt
 */
#pragma once

#include <functional>
#include <unordered_map>

namespace HomologyInference
{

template <class T> using HashFunction = std::function<size_t(const T&)>;

template <class T>
inline void hash_combine(std::size_t& hash, const T& v)
{
    // https://stackoverflow.com/questions/2590677/how-do-i-combine-hash-values-in-c0x
    std::hash<T> hasher;
    hash ^= hasher(v) + 0x9e3779b9 + (hash<<6) + (hash>>2);
}

template <class RangeT>
inline size_t range_hash(const RangeT _range)
{
    size_t hash = 0;

    for (const auto& e : _range)
        hash_combine(hash, e);

    return hash;
}

struct pair_hash
{
    template <class T1, class T2>
    std::size_t operator() (const std::pair<T1, T2>& pair) const
    {
        size_t hash = 0;
        hash_combine(hash, pair.first);
        hash_combine(hash, pair.second);
        return hash;
    }
};

}
