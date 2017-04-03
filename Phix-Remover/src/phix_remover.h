/*The idea of this program is to create a quick kmer lookup table for phix.
 * It accomplishes this by making 2 arrays of 1<<16 size, because we represent
 * each possible kmer in the 2 bit format (A -> 00, T -> 11, C -> 01, G->10) and
 * need foward and Reverse strands. Each read is then parsed into 8mers (again
 * represented by their 2 bit format, and checked against these lookup tables
 * (simply the arrays)*/ 

#ifndef AT_TRIM_H
#define AT_TRIM_H

#define BOOST_DYNAMIC_BITSET_DONT_USE_FRIENDS

#include "ioHandler.h"
#include <map>
#include <unordered_map>
#include <boost/dynamic_bitset.hpp>
#include <boost/functional/hash.hpp>
#include <algorithm>
#include <bitset>
#include "utils.h"
#include <unordered_set>

/*This will be the size_t will be the number of hits to that kmer
 * and the number of entries will be the number of possible kmers
 * ~64k. This will create a quick lookup table to use for both the
 * forward and reverse complement of phix.*/
//typedef std::array< size_t, 1<<(kmerBits) > kmerArray;

class dbhash {
public:
    std::size_t operator() ( const boost::dynamic_bitset<>& bs) const {
        return boost::hash_value(bs.m_bits);
    }
};

typedef std::unordered_set < boost::dynamic_bitset<>, dbhash> kmerSet;

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


int setter(boost::dynamic_bitset<> &lookup, size_t loc, char c, bool rc) {

    if (c == 'A') {
        lookup[loc] = (0 ^ rc);
        lookup[loc+1] = (0 ^ rc);
    } else if (c == 'T') {
        lookup[loc] = (1 ^ rc);
        lookup[loc+1] = (1 ^ rc);
    } else if (c == 'C') {
        lookup[loc] = (1 ^ rc);
        lookup[loc+1] = (0 ^ rc);
    } else if (c == 'G') {
        lookup[loc] = (0 ^ rc)  ;
        lookup[loc+1] = (1 ^ rc);
    } else {
        return 0; //N
    }
    return 1;
}

/*Check read will only return back the number of hits to the lookup tables,
 * the look up tables at this point should all ready be set for phix at this point. It 
 * will return the maximum number of hits from either the forward or rc 
 * lookup table. If there are enough hits, then the read will be assumed to be phix,
 * and will be discarded.*/
size_t check_read(const Read &r, const kmerSet &lookup, const kmerSet &lookup_rc, size_t kmerSize) {
    
    boost::dynamic_bitset <> forwardBits(kmerSize * 2);
    size_t kmerCheck = kmerSize - 1;
    size_t bin = 0;
    size_t hits = 0, hits_rc = 0;
    size_t loc = (kmerSize * 2) - 2;
    std::string seq = r.get_seq();
    
    for (std::string::iterator bp = seq.begin(); bp < seq.end(); ++bp) {

        if ( !setter(forwardBits, loc, *bp, false) ) {
            kmerCheck = kmerSize - 1;
            forwardBits.reset();
            loc = (kmerSize * 2) -2 ;
        } else {
            if (!kmerCheck) {
                hits += lookup.find(forwardBits) != lookup.end();
                hits_rc += lookup_rc.find(forwardBits) != lookup_rc.end();
                forwardBits <<= 2;
            } else {
                --kmerCheck;
                loc -= 2;
            }
        }
    } 
    
    return std::max(hits, hits_rc);
}


template <class T, class Impl>
void helper_discard(InputReader<T, Impl> &reader, std::shared_ptr<OutputWriter> pe, std::shared_ptr<OutputWriter> se, Counter& c, kmerSet &lookup, kmerSet &lookup_rc, size_t hits, bool checkR2, size_t kmerSize) {
    
    while(reader.has_next()) {
        auto i = reader.next();
        PairedEndRead* per = dynamic_cast<PairedEndRead*>(i.get());        
        if (per) {
            size_t val = check_read(per->get_read_one(), lookup, lookup_rc, kmerSize);
            if (checkR2) {
                size_t val2 = check_read(per->get_read_two(), lookup, lookup_rc, kmerSize);
                val = std::max(val, val2);
            }

            if (val < hits) {
                writer_helper(per, pe, se, false, c);
            } else {
                ++c["R1_Discarded"];
                ++c["R2_Discarded"];
            }

        } else {
            SingleEndRead* ser = dynamic_cast<SingleEndRead*>(i.get());
            
            if (ser) {
                size_t val = check_read(ser->get_read(), lookup, lookup_rc, kmerSize);
                if (val < hits) {
                    writer_helper(ser, pe, se, false, c);
                } else {
                    ++c["SE_Discarded"];
                }
            } else {
                throw std::runtime_error("Unknown read type");
            }
        }
    }
}

/*Hand phix read to set lookup tables
 *The lookup tables might make more sense in a map format
 * the key would be the two bit format of the 8mer (so a uint_16_t)
 * and the value would be number of hits
 * The reason we are using arrays though is because we can use the
 * pigeon hole theorm and have a constant lookup time because
 * the arrays are small enough.
 * These lookup tables give us a quick way to count the "hits" of 8mers
 * from the sequencing read. If there are enough "hits", we can assume it
 * is a phix read.
 * */
void setLookup(kmerSet &lookup, kmerSet &lookup_rc, Read &rb, size_t kmerSize) {
    size_t kmerCheck = kmerSize - 1;
    size_t bitKmer = kmerSize * 2;
    boost::dynamic_bitset <> forwardBits(kmerSize * 2);
    boost::dynamic_bitset <> reverseBits(kmerSize * 2);
    size_t loc = (kmerSize * 2) - 2;
    size_t loc_rc = 0;
    std::string seq = rb.get_seq();

    for (std::string::iterator bp = seq.begin(); bp < seq.end(); ++bp) {
        /*Cannot handle Ns in 2 bit form
         * reset kmerCheck and bits*/

        if ( !setter(forwardBits, loc, *bp, false) ) {
            kmerCheck = kmerSize - 1;
            forwardBits.reset();
            reverseBits.reset();
            loc = (kmerSize * 2) - 2;
            loc_rc = 0;
        } else {
            setter( reverseBits, loc_rc, *bp, true);
            if (!kmerCheck) {
                lookup.insert(forwardBits);
                lookup_rc.insert(reverseBits);
                //++lookup[forwardBits.to_ulong()];
                forwardBits <<= 2;
                reverseBits >>= 2;
                //++lookup_rc[reverseBits.to_ulong()];
            } else {
                --kmerCheck;
                loc -= 2;
                loc_rc += 2;
            }
        }
        
        /* This is for the reverse complement of the lookup string.
         * It is possible to either see the forward or RC strand of Phix
         * Versus checking every forward and reverse strand of each sequence seen
         * it is quicker to create an RC lookup table and check both the forward and cehck both.
         * Shifts the 0 and 1st bit off. Then use the binary rc of the bp to set the 14 and 15 bit*/

        /*If kmerCheck is zero (all bp accounted for), add a point to the lookup table*/
    }
    
}

#endif

