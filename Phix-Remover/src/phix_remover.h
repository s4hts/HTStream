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
#include <bitset>

typedef std::unordered_map<std::string, size_t> Counter;
const size_t kmer = 8;
const size_t kmerBits = kmer*2;
typedef std::array< size_t, 1<<(2*kmer) > kmerArray;

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

size_t check_read(const Read &r, kmerArray lookup, kmerArray lookup_rc) {
    
    std::bitset <kmerBits> forwardBits;
    size_t kmerCheck = kmer - 1;
    size_t bin = 0;
    size_t hits = 0, hits_rc = 0;

    std::string seq = r.get_seq();
    
    for (std::string::iterator bp = seq.begin(); bp < seq.end(); ++bp) {
        
        if ((bin = getBin(*bp)) == 5) { //bp = N
            kmerCheck = kmer - 1;
            forwardBits.reset();
        }

        forwardBits &= ~(3 << 14); //only check forward strand against both RC and non-rc
        forwardBits <<= 2;
        forwardBits ^= bin;
        
        if (!kmerCheck) {
            if (lookup[forwardBits.to_ulong()]) {
                ++hits;
            } else if (lookup_rc[forwardBits.to_ulong()]) {
                ++hits_rc;
            }
        } else {
            --kmerCheck;    
        }
    } 
    return std::max(hits, hits_rc);
}


template <class T, class Impl>
void helper_discard(InputReader<T, Impl> &reader, std::shared_ptr<OutputWriter> pe, std::shared_ptr<OutputWriter> se, Counter& counters, kmerArray &lookup, kmerArray &lookup_rc, size_t hits, bool checkR2) {
    
    while(reader.has_next()) {
        auto i = reader.next();
        ++counters["TotalRecords"];
        PairedEndRead* per = dynamic_cast<PairedEndRead*>(i.get());        
        if (per) {
            size_t val = check_read(per->get_read_one(), lookup, lookup_rc);
            if (checkR2) {
                size_t val2 = check_read(per->get_read_two(), lookup, lookup_rc);
                val = std::max(val, val2);
            }

            if (val < hits) {
                pe->write(*per);
            }

        } else {
            SingleEndRead* ser = dynamic_cast<SingleEndRead*>(i.get());
            
            if (ser) {
                size_t val = check_read(ser->get_read(), lookup, lookup_rc);
                if (val > hits) {
                    se->write(*ser);
                }
            } else {
                throw std::runtime_error("Unknown read type");
            }
        }
    }
}

/*Hand phix read to set lookup tables*/
void setLookup(kmerArray &lookup, kmerArray &lookup_rc, Read &rb) {
    size_t kmerCheck = kmer - 1;
    size_t bin = 0;
    size_t bin_rc = 0;

    std::bitset <kmerBits> forwardBits;
    std::bitset <kmerBits> reverseBits;

    std::string seq = rb.get_seq();

    for (int i = 0; i < 1<<(kmerBits); ++i) {
        lookup[i] = 0;
        lookup_rc[i] = 0;
    }
    
    for (std::string::iterator bp = seq.begin(); bp < seq.end(); ++bp) {
        /*Cannot handle Ns in 2 bit form
         * reset kmerCheck and bits*/
        if ((bin = getBin(*bp)) == 5) { 
            kmerCheck = kmer - 1;
            forwardBits.reset();
            reverseBits.reset();
        }
        /*Takes the first two bits and nots them (RC the bp)*/
        bin_rc = (bin ^ ((1 << 2) - 1));

        /*Sets the 14 and 15th bit to zero to allow room to shift
         * over by two. The 0 and 1st bit are not able to be set
         * based on the bp seen*/
        forwardBits &= ~(3 << ((kmerBits)-2));
        forwardBits <<= 2;
        forwardBits ^= bin;
        
        /* This is for the reverse complement of the lookup string.
         * It is possible to either see the forward or RC strand of Phix
         * Versus checking every forward and reverse strand of each sequence seen
         * it is quicker to create an RC lookup table and check both the forward and cehck both.
         * Shifts the 0 and 1st bit off. Then use the binary rc of the bp to set the 14 and 15 bit*/
        reverseBits >>=  2;
        reverseBits ^= (bin_rc << 14);

        /*If kmerCheck is zero (all bp accounted for), add a point to the lookup table*/
        if (!kmerCheck) {
            ++lookup[forwardBits.to_ulong()];
            ++lookup_rc[reverseBits.to_ulong()];
        } else {
            --kmerCheck;
        }
    }
    
}

#endif

