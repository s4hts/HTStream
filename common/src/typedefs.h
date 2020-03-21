#ifndef TYPDEF_H
#define TYPDEF_H

#include <cstdint>
#include <map>

template <typename G>
using Generic = std::tuple <const std::string, G>;

typedef std::tuple <const std::string, uint64_t&> Label;
typedef std::tuple <const std::string, std::string> sLabel;
typedef std::tuple<uint_fast64_t, uint_fast64_t> Vector;

typedef std::vector<uint_fast64_t> Vec;
typedef std::vector<Vec> Mat;

#endif
