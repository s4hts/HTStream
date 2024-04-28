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

const uint64_t QUAL_MAX = 94;
const uint64_t DEFAULT_QUAL_OFFSET = 33;
const std::vector<char> READ_OPTIONS{'F', 'R', 'B'}; // possible read options
const std::vector<char> DEL_OPTIONS{'-', '_', ':', '\0'}; // possible char options (last one is to allow unset)


#endif
