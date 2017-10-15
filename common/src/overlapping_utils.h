#ifndef OVERLAPPING_UTILS_H
#define OVERLAPPING_UTILS_H

#include "ioHandler.h"
#include <map>
#include <unordered_map>
#include <boost/dynamic_bitset.hpp>
#include <boost/functional/hash.hpp>
#include <algorithm>
#include <bitset>
#include <utility>
#include "utils.h"

typedef std::unordered_multimap<std::string, std::size_t> seqLookup;

unsigned int getInsertLength(PairedEndRead &pe , const double misDensity, const size_t &minOver, const size_t &checkLengths, const size_t kmer, const size_t kmerOffset, const size_t minLength) ;

seqLookup readOneMap(std::string seq1, const size_t kmer, const size_t kmerOffset) ;

unsigned int getOverlappedReads(Read &r1, Read &r2, const seqLookup &seq1Map,  const double misDensity, const size_t &minOver, const size_t &checkLengths, const size_t kmer);

//previously check_read
unsigned int getInsertSize(Read &r1, Read &r2, size_t loc1, size_t loc2, const double misDensity, size_t minOverlap);

unsigned int checkIfOverlap(Read &r1, Read &r2, size_t loc1, size_t loc2, const double misDensity, size_t minOverlap);


unsigned int getInsertSize(PairedEndRead &pe , const double misDensity, const size_t &minOver, const size_t &checkLengths, const size_t kmer, const size_t kmerOffset, const size_t minLength);

#endif
