#ifndef TYPDEF_H
#define TYPDEF_H

#include <cstdint>
#include <map>

typedef std::tuple <const std::string, uint64_t&> Label;
typedef std::tuple<uint_fast64_t, uint_fast64_t> Vector;

// seqMatrix is for stats counters, both seq character and seq quality
// pos, character, count
typedef std::tuple<uint_fast64_t, const char, uint_fast64_t>  ;

#endif
