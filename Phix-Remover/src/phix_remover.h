#ifndef AT_TRIM_H
#define AT_TRIM_H
//  this is so we can implment hash function for dynamic_bitset
#define BOOST_DYNAMIC_BITSET_DONT_USE_FRIENDS

#include "ioHandler.h"
#include <map>
#include <unordered_map>
#include <boost/dynamic_bitset.hpp>
#include <boost/functional/hash.hpp>
#include <algorithm>

typedef std::unordered_map<std::string, size_t> Counter;

uint8_t getBin(char c) {
    if (c == 'A') 
        return 0;
    else if (c == 'T')
        return 3;
    else if (c == 'C')
        return 1;
    else if (c == 'G')
        return 2;
    return 5;
}

size_t check_read(const Read &r, std::array<size_t, 1<<(2*8) > lookup, std::array< size_t, 1 << (2*8) > lookup_rc, size_t kmerSize) {
    size_t index = 0;
    size_t kmer = kmerSize - 1;
    size_t bin = 0;
    size_t hits = 0, hits_rc = 0;
    std::string seq = r.get_seq();
    
    for (std::string::iterator bp = seq.begin(); bp < seq.end(); ++bp) {
        
        if ((bin = getBin(*bp)) == 5) { //bp = N
            kmer = kmerSize - 1;
            index = 0;
        }

        index &= ~(3 << 14);
        index <<= 2;
        index ^= bin;
        
        if (!kmer) {
            if (lookup[index]) {
                ++hits;
            } else if (lookup_rc[index]) {
                ++hits_rc;
            }
        } else {
            --kmer;    
        }
    } 
    return std::max(hits, hits_rc);
}


template <class T, class Impl>
void helper_discard(InputReader<T, Impl> &reader, std::shared_ptr<OutputWriter> pe, std::shared_ptr<OutputWriter> se, Counter& counters, std::array<size_t, 1<<(2*8) > lookup, std::array<size_t, 1<<(2 * 8)> lookup_rc, size_t kmerSize, size_t hits, bool checkR2) {
    
    while(reader.has_next()) {
        auto i = reader.next();
        ++counters["TotalRecords"];
        PairedEndRead* per = dynamic_cast<PairedEndRead*>(i.get());        

        if (per) {
            size_t val = check_read(per->get_read_one(), lookup, lookup_rc, kmerSize);
            if (checkR2) {
                size_t val2 = check_read(per->get_read_two(), lookup, lookup_rc, kmerSize);
                val = std::min(val, val2);
            }

            if (val > hits) {
                pe->write(*per);
            }

        } else {
            SingleEndRead* ser = dynamic_cast<SingleEndRead*>(i.get());
            
            if (ser) {
                size_t val = check_read(ser->get_read(), lookup, lookup_rc, kmerSize);
                if (val > hits) {
                    se->write(*ser);
                }
            } else {
                throw std::runtime_error("Unknown read type");
            }
        }
    }

}

void setLookup(std::array<size_t, 1<<2*8> &lookup, std::array<size_t, 1<<2*8> &lookup_rc, Read &rb, size_t kmerSize) {
    size_t index = 0;
    size_t index_rc = 0;
    size_t kmer = kmerSize - 1;
    size_t bin = 0;
    size_t bin_rc = 0;

    std::string seq = rb.get_seq();

    for (int i = 0; i < 1<<(2*kmerSize); ++i) {
        lookup[i] = 0;
        lookup_rc[i] = 0;
    }
    
    for (std::string::iterator bp = seq.begin(); bp < seq.end(); ++bp) {
        if ((bin = getBin(*bp)) == 5) { //bp = N
            kmer = kmerSize - 1;
            index = 0;
            index_rc = 0;
        }
        bin_rc = (bin ^ ((1 << 2) - 1));
        index &= ~(3 << 14);
        index <<= 2;
        index ^= bin;
        
        index_rc >>=  2;
        index_rc ^= (bin_rc << 14);
        if (!kmer) {
            ++lookup[index];
            ++lookup_rc[index_rc];
            std::cout << index_rc << '\n';
        } else {
            --kmer;    
        }
    }
    
}

#endif

